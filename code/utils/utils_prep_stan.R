## ----------------------------------------------------------------
##                  Stan data list preparation                   -
## ----------------------------------------------------------------

compute_backward_delay <- function(time_series, forward_delay_dist){
  bd <- matrix(0, nrow = length(time_series) + length(forward_delay_dist) - 1, ncol = length(forward_delay_dist))
  for (i in 1:length(time_series)){
    for (d in 1:length(forward_delay_dist)){
      bd[i + d - 1, d] <- bd[i + d - 1, d] + time_series[i] * forward_delay_dist[d]
    }
  }
  for (i in 1:(length(time_series) + length(forward_delay_dist) - 1)){
    if (sum(bd[i, ])>0) bd[i, ] <- bd[i, ] / sum(bd[i, ])
  }
  return(bd)
}

random_inits <- function(stan_data_list, model_type) {
  return (NULL)
}

default_inits <- function(stan_data_list, model_type) {
  known_by_reference_date <- rowSums(stan_data_list$reported_known)
  if (sum(stan_data_list$reported_unknown)>0) {
    unknown_by_reference_date <- rep(0, stan_data_list$T)
    for (i in 1:stan_data_list$T) {
      # assume flat reporting delay
      unknown_fraction <- stan_data_list$reported_unknown[i] / (stan_data_list$D + 1)
      unknown_by_reference_date[max((i - stan_data_list$D), 1):i] <- unknown_by_reference_date[max((i - stan_data_list$D), 1):i] + unknown_fraction
    }
    all_by_reference_date <- unknown_by_reference_date + known_by_reference_date
    alpha_emp <- known_by_reference_date / all_by_reference_date
    alpha_emp <- mean(alpha_emp[1:7])
  } else {
    all_by_reference_date <- known_by_reference_date
    alpha_emp <- rep(1-1e-4,7)
  }
  alpha_emp_logit <- qlogis(alpha_emp)

  # get crude estimate of baseline hazard, assume beta and eta are zero
  gamma_emp <- rowMeans(sapply(1:7, function(i) stan_data_list$reported_known[i, ] / sum(stan_data_list$reported_known[i, ])))
  gamma_emp_haz <- get_hazard_from_p(gamma_emp)
  gamma_emp_haz_reduced <- rep(0, stan_data_list$n_delays)
  for (i in 1:stan_data_list$n_delays) {
    gamma_emp_haz_reduced[i] <- mean(gamma_emp_haz[stan_data_list$delay_idx == i])
  }

  if (!(model_type %in% c("base","base_old"))) {
  # assume infections = symptom onsets shifted by mean incubation period
  medianInc <- weighted.median(0:(length(stan_data_list$latent_delay_dist) - 1), stan_data_list$latent_delay_dist)
  shift_backward <- function(x, shift) c(x[(1 + shift):length(x)], rep(NA, times = shift))
  I_length <- stan_data_list$L + stan_data_list$n_lambda_pre + stan_data_list$T
  start_t <- stan_data_list$L + stan_data_list$n_lambda_pre
  I_by_reference_date <- rep(NA, I_length)
  I_by_reference_date[(start_t + 1):(start_t + length(all_by_reference_date))] <- all_by_reference_date
  firstI <- I_by_reference_date[start_t + 1]
  lastI <- I_by_reference_date[start_t + length(all_by_reference_date)]
  I_by_reference_date <- shift_backward(I_by_reference_date, medianInc)
  I_by_reference_date[1:(start_t - medianInc)] <- firstI
  I_by_reference_date[(start_t + length(all_by_reference_date) - medianInc):(start_t + length(all_by_reference_date))] <- lastI
  # back-calculate iota values and R time series
  gt_mean <- weighted.mean(0:(length(stan_data_list$generation_time_dist) - 1), stan_data_list$generation_time_dist)
  gt_var <- weighted.mean((0:(length(stan_data_list$generation_time_dist) - 1))^2, stan_data_list$generation_time_dist) - gt_mean^2
  R_emp <- c(NA, NA, EpiEstim::estimate_R(
    incid = I_by_reference_date,
    method = "parametric_si",
    config = EpiEstim::make_config(
      list(
        mean_si = gt_mean,
        std_si = sqrt(gt_var),
        t_start = 2:(length(I_by_reference_date) - 1),
        t_end = 3:length(I_by_reference_date),
        mean_prior = 1
      )
    )
  )$R$`Mean(R)`)
  }

  return(
    function() {
      inits <- list()

      # use an empirical estimate of alpha and assume it remains roughly unchanged over time
      inits[["alpha_logit_start"]] <- alpha_emp_logit
      inits[["alpha_logit_noise"]] <- rep(0, stan_data_list$T)
      inits[["alpha_logit_sd"]] <- rtnorm(n = 1, mean = stan_data_list$alpha_logit_sd_prior_mu, sd = stan_data_list$alpha_logit_sd_prior_sd / 4, a = 0) # truncated

      inits[["gamma"]] <- gamma_emp_haz_reduced
      inits[["beta"]] <- rep(1e-4, stan_data_list$n_beta)
      inits[["eta"]] <- rep(1e-4, stan_data_list$n_eta)

      inits[["xi_negbinom"]] <- rtnorm(n = 1, mean = stan_data_list$xi_negbinom_prior_mu, sd = stan_data_list$xi_negbinom_prior_sd / 4, a = 0)

      if (model_type %in% c("base","base_old")) {
        inits[["ets_alpha"]] <- 1 - 1e-4
        inits[["ets_beta"]] <- 1 - 1e-4
        inits[["ets_phi"]] <- 1 - 1e-4
        inits[["lambda_log_start_values"]] <- c(log(stan_data_list$expected_cases_start),rep(0,1+stan_data_list$ets_diff))
        inits[["lambda_log_noise"]] <- rep(0, stan_data_list$T + stan_data_list$n_lambda_pre - 1 - stan_data_list$ets_diff) #rnorm(n = stan_data_list$T + stan_data_list$n_lambda_pre - 1)
        inits[["lambda_log_sd"]] <- stan_data_list$lambda_log_sd_prior_mu + stan_data_list$lambda_log_sd_prior_sd / 2 #rtnorm(n = 1, stan_data_list$lambda_log_sd_prior_mu, sd = stan_data_list$lambda_log_sd_prior_sd / 4, a = 0) # truncated
      } else if (model_type == "latent") {
        # use initial expected cases as proxy for initial expected infections
        inits[["iota_log_start"]] <- rnorm(n = 1, mean = log(stan_data_list$expected_cases_start), sd = 0.4)
        inits[["iota_log_sd"]] <- 1
        inits[["iota_log_noise"]] <- c(0, diff(log(I_by_reference_date))[1:(stan_data_list$T + stan_data_list$L + stan_data_list$n_lambda_pre)])
      } else if (model_type == "renewal") {
        inits[["ets_alpha"]] <- 1 - 1e-4 # rbeta(n = 1, stan_data_list$ets_alpha_prior_alpha, stan_data_list$ets_alpha_prior_beta)
        inits[["ets_beta"]] <- 1 - 1e-4 # rbeta(n = 1, stan_data_list$ets_beta_prior_alpha, stan_data_list$ets_beta_prior_beta)
        inits[["ets_phi"]] <- 1 - 1e-4 # rbeta(n = 1, stan_data_list$ets_phi_prior_alpha, stan_data_list$ets_phi_prior_beta)
        inits[["R_level_start"]] <- R_emp[stan_data_list$max_gen + 1] # rtnorm(n = 1, mean = stan_data_list$R_level_start_prior_mu, sd = stan_data_list$R_level_start_prior_sd / 4, a = 0) # truncated
        inits[["R_trend_start"]] <- 1e-4 # rnorm(n = 1, mean = stan_data_list$R_trend_start_prior_mu, sd = stan_data_list$R_trend_start_prior_sd / 4)
        inits[["R_sd"]] <- 1
        inits[["R_noise"]] <- rep(0, stan_data_list$L + stan_data_list$n_lambda_pre + stan_data_list$T - stan_data_list$max_gen - 1)# c(0, diff(R_emp))[(stan_data_list$max_gen + 2):(stan_data_list$L + stan_data_list$n_lambda_pre + stan_data_list$T)]
        # use initial expected cases as proxy for infections / expected infections
        inits[["iota_log_ar_start"]] <- log(I_by_reference_date[1])
        inits[["iota_log_ar_sd"]] <- 1
        inits[["iota_log_ar_noise"]] <- c(0, diff(log(I_by_reference_date))[1:(stan_data_list$max_gen - 1)])
        inits[["I"]] <- I_by_reference_date
      }
      return(inits)
    }
  )
}

prepare_data_list <- function(prep_data_complete,
                              prep_data_missing,
                              now,
                              start_date,
                              holidays = NULL,
                              D,
                              delay_idx,
                              n_lambda_pre,
                              model_type = "base") {

  # --------------------------------------------
  # Reporting triangle
  reporting_matrices <- prepare_reporting_data(prep_data_complete, prep_data_missing, now, start_date, D)

  # estimate initial expected cases based on reported known + unknown assuming a flat delay distribution
  expected_cases_start <- sum(reporting_matrices[["reported_known"]][1, ]) +
    mean(reporting_matrices[["reported_unknown"]][1:D])

  # ---------------------------
  # Covariates for reporting delay
  all_dates <- seq(start_date, now, by = "1 day")
  tswp <- seq(start_date - n_lambda_pre, now, by = "1 day") # time series including the days of the modeling phase
  T <- length(all_dates) - 1

  ## 1) Segmented linear model for change over time
  ddChangepoint <- rev(seq(now, start_date, by = "-2 weeks"))
  Z_idx <- rev(rep(length(ddChangepoint):1, each = 14)[1:(1 + T)])
  Z <- array(NA, dim = c(T + 1 + n_lambda_pre, length(ddChangepoint)), dimnames = list(as.character(tswp), 1:length(ddChangepoint)))
  for (t_pre in 1:n_lambda_pre) {
    Z[t_pre, ] <- rep(0, length(ddChangepoint)) # we don't have a changepoint model before the inference phase, but fill the space up with 0s to have the same dimensionality as the other matrices
  }
  for (t in 0:T) {
    covariates_ind_days <- c(0, as.numeric(all_dates <= (all_dates[t + 1]))[-1])
    Z[n_lambda_pre + 1 + t, ] <- tapply(covariates_ind_days, Z_idx, sum) # add up counts in each changepoint interval
  }
  # scale columns by maximum for better efficiency
  Z <- Z / (max(Z) * 0.68)


  ## 2) Weekday effect
  # Make designmatrix for weekdays
  # Make a factor for each day of the week
  wdays <- levels(lubridate::wday(tswp, label = TRUE))[-1]
  # Create the extension of the W matrix
  W <- array(NA,
    dim = c(length(tswp), D, length(wdays)),
    dimnames = list(as.character(tswp), paste("delay", 0:(D - 1), sep = ""), wdays)
  )

  # Replace holidays with sunday indicator
  get_wdays <- function(nc_dates_t, D) {
    wdays <- lubridate::wday(nc_dates_t + 0:(D - 1), label = TRUE)
    if (!is.null(holidays)) {
      wdays[(nc_dates_t + 0:(D - 1)) %in% holidays] <- "Sun"
    }
    return(wdays)
  }

  # Loop over all times and lags
  for (t in seq_len(length(tswp))) {
    for (w in seq_len(length(wdays))) {
      #
      W[t, , w] <- as.numeric(get_wdays(tswp[t], D) == wdays[w])
    }
  }

  # ---------------------------
  # Incubation period distribution
  L <- 14
  latent_delay_dist <- get_incubation_dist(gamma_mean = 5.3, gamma_sd = 3.2, maxInc = L)

  # ---------------------------
  # Generation time distribution
  maxGen <- 10
  generation_time_dist <- get_generation_dist(gamma_mean = 4.8, gamma_sd = 2.3, maxGen = maxGen)

  # ---------------------------
  # Reduce design matrices
  if (!all(delay_idx == seq(1, D + 2))) {
    W <- reduce_design_matrix(W, delay_idx)
  }

  #-------------------------------------------
  # Combine all data/preprocessed objects into list

  stan_data_list <- list() # directly used in the Stan model
  additional_info <- list() # used by other parts of the program or as info for the user

  # durations
  stan_data_list[["T"]] <- T + 1
  stan_data_list[["D"]] <- D
  stan_data_list[["Z"]] <- Z
  stan_data_list[["W"]] <- W
  stan_data_list[["n_lambda_pre"]] <- n_lambda_pre
  pre_diff <- D - n_lambda_pre

  # observed data
  stan_data_list[["reported_known"]] <- reporting_matrices[["reported_known"]]
  stan_data_list[["reported_unknown"]] <- reporting_matrices[["reported_unknown"]]
  stan_data_list[["expected_cases_start"]] <- expected_cases_start
  if (pre_diff>0) stan_data_list[["emp_alpha"]] <- reporting_matrices[["emp_alpha"]][1:pre_diff]

  # reporting model
  stan_data_list[["n_delays"]] <- max(delay_idx) - 2
  stan_data_list[["delay_idx"]] <- delay_idx
  stan_data_list[["n_beta"]] <- length(ddChangepoint)
  stan_data_list[["n_eta"]] <- length(wdays)

  # latent process details
  if (model_type == "latent" | model_type %in% c("renewal", "renewal_noncentered", "renewal_deterministic")) {
    stan_data_list[["L"]] <- L
    stan_data_list[["latent_delay_dist"]] <- latent_delay_dist
  }

  # renewal process details
  if (model_type %in% c("renewal", "renewal_noncentered", "renewal_deterministic")) {
    stan_data_list[["max_gen"]] <- maxGen
    stan_data_list[["generation_time_dist"]] <- generation_time_dist
  }

  # all dates
  additional_info[["all_dates"]] <- all_dates

  return(list(stan_data_list = stan_data_list, additional_info = additional_info))
}

## ----------------------------------------------------------------
##              Preparation of reporting triangles               -
## ----------------------------------------------------------------

prepare_reporting_data <- function(prep_data_complete, prep_data_missing, now, start_date, D) {
  if (!("n" %in% names(prep_data_complete))) prep_data_complete %<>% mutate(n = 1)
  if (!("n" %in% names(prep_data_missing))) prep_data_missing %<>% mutate(n = 1)

  # -------------------------
  # Cases with observed type 1 event
  reporting_triangle <- prep_data_complete %>%
    filter(event1_date >= start_date, event2_date <= now) %>%
    mutate(delay = as.numeric(event2_date - event1_date))

  # Currently, we exclude observations with negative delay. If alternative behavior is wanted (e.g. trim to zero, or set missing),
  # this should already be done before calling this function. Observations with delay longer than the specified maximum delay are truncated.
  n_negative_delay <- reporting_triangle %>%
    ungroup() %>%
    filter(delay < 0) %>%
    summarize(n = sum(n), .groups = "drop") %>%
    pull(n)
  n_maximum_delay <- reporting_triangle %>%
    ungroup() %>%
    filter(delay > D) %>%
    summarize(n = sum(n), .groups = "drop") %>%
    pull(n)
  if (n_negative_delay > 0) warning(paste("Removed ", n_negative_delay, " case(s) with negative delay.", sep = ""))
  if (n_maximum_delay > 0) warning(paste("Truncated ", n_maximum_delay, " cases with delay > ", D, ".", sep = ""))

  reporting_triangle %<>%
    filter(delay >= 0) %>% # remove negative delays
    mutate(delay = ifelse(delay > D, D, delay)) # truncate delays longer than D

  # Compute triangle
  reporting_triangle %<>%
    group_by(event1_date, delay) %>%
    dplyr::summarize(n = sum(n), .groups = "drop") %>%
    complete(delay = 0:D, event1_date, fill = list(n = 0)) %>%
    arrange(delay, event1_date) %>%
    pivot_wider(names_from = delay, values_from = n) %>%
    ungroup() %>%
    complete(event1_date = seq.Date(start_date, now, by = "1 day")) %>%
    mutate(across(-event1_date, replace_na, replace = 0))

  # Convert to matrix
  reporting_matrix <- as.matrix(reporting_triangle %>% select(-event1_date))
  dimnames(reporting_matrix)[[1]] <- as.character(reporting_triangle$event1_date)

  # -------------------------
  # Cases with unobserved type 1 event
  unobserved <- prep_data_missing %>%
    group_by(event2_date) %>%
    dplyr::summarize(n = sum(n), .groups = "drop") %>%
    filter(event2_date >= start_date, event2_date <= now) %>%
    complete(event2_date = seq.Date(start_date, now, by = "1 day"), fill = list(n = 0))

  unobserved_matrix <- unobserved %>% pull(n)

  # -------------------------
  # Seed for expected cases with unobserved type 1 event (based on backward delays)
  emp_backward_delay_dist <- prep_data_complete %>%
    filter(event2_date >= start_date, event2_date <= now) %>%
    get_empirical_backward_delay(D)

  emp_alpha <- unobserved %>%
    inner_join(emp_backward_delay_dist, by = "event2_date") %>%
    mutate(a = n * p, event1_date = event2_date - delay) %>%
    group_by(event1_date) %>%
    summarize(unobserved = sum(a)) %>%
    inner_join(reporting_triangle %>%
      transmute(event1_date, observed = rowSums(across(-event1_date))),
    by = "event1_date"
    ) %>%
    mutate(alpha = observed / (observed + unobserved)) %>%
    ungroup() %>% 
    mutate(alpha = ifelse(is.nan(alpha),NA,alpha)) %>% 
    fill(alpha, .direction = "downup") %>% 
    pull(alpha)

  return(list(reported_known = reporting_matrix, reported_unknown = unobserved_matrix, emp_alpha = emp_alpha))
}

get_empirical_backward_delay <- function(df, D) {
  df %>% 
    mutate(delay = as.numeric(event2_date - event1_date)) %>%
    filter(delay >= 0) %>% # remove negative delays
    mutate(delay = ifelse(delay > D, D, delay)) %>% # truncate delays longer than D
    count(event2_date, delay) %>%
    complete(event2_date, delay = 0:D, fill = list(n = 0)) %>%
    group_by(event2_date) %>%
    transmute(delay, p = n / sum(n)) %>% 
    return()
}

## -------------------------------------------------------------
##  Reduction of design matrices for lower resolution intervals
## -------------------------------------------------------------

reduce_design_matrix <- function(design_matrix, delay_idx) {
  n_time_points <- dim(design_matrix)[1]
  n_reduced_times <- max(delay_idx) - 2
  n_covariates <- dim(design_matrix)[3]
  reduced_matrix <- array(dim = c(n_time_points, n_reduced_times, n_covariates))

  for (d in 1:n_reduced_times) {
    delay_indices <- delay_idx[1:(length(delay_idx) - 2)] == d
    if (sum(delay_indices) == 1) {
      reduced_matrix[, d, ] <- design_matrix[, delay_indices, ]
    } else {
      reduced_matrix[, d, ] <- t(sapply(1:n_time_points, function(t) colMeans(design_matrix[t, delay_indices, ])))
    }
  }
  return(reduced_matrix)
}

## ----------------------------------------------------------------
##                Prior for the delay distribution               -
## ----------------------------------------------------------------

## ------------------------------------
##  Dirichlet gamma prior
## ------------------------------------

# The goal of the below functions is to find a good prior (through mu and sigma of a normal distribution)
# for the gamma parameter of the model (baseline logit-hazard for delay d).
# This is achieved by choosing the mu and sigma such that the resulting delay distribution
# is as similar as possible to a Dirichlet distribution where all days within the horizon D have the same probability.
# The kappa parameter (related to the variance of the Dirichlet) can be chosen by the user.

# Code taken from Günther et. al (2021), see https://github.com/FelixGuenther/nc_covid19_bavaria,
# with minor improvement to respect non-negativity of variance parameter

expectation_p_haz_plogisN <- function(mu, sigma, p_smaller) {
  f <- function(x) {
    (1 - p_smaller) * plogis(x) * dnorm(x, mu, sigma)
  }
  int <- integrate(f, lower = -Inf, upper = Inf)
  int$value
}

# Variance of P(D>=d)*P(D=d|D>=d), with P(D=d|D>=d)=1/(1+exp(-x)) and X~N(mu, sigma^2)
var_p_haz_plogisN <- function(mu, sigma, p_smaller) {
  E <- expectation_p_haz_plogisN(mu, sigma, p_smaller)
  f <- function(x) {
    (((1 - p_smaller) * plogis(x)) - E)^2 * dnorm(x, mu, sigma)
  }
  int <- integrate(f, lower = -Inf, upper = Inf)
  int$value
}

get_prior_gamma <- function(gd.prior.kappa = 1, D = 20) {
  alpha_i <- gd.prior.kappa
  alpha_0 <- (D + 1) * gd.prior.kappa
  e_dir <- alpha_i / alpha_0
  var_dir <- (alpha_i * (alpha_0 - alpha_i)) / (alpha_0^2 * (alpha_0 + 1))

  sum_sq_diff <- function(theta, p_smaller) {
    mu <- theta[1]
    sigma <- exp(theta[2]) + 10^-10
    E_diff <- expectation_p_haz_plogisN(mu, sigma, p_smaller) - e_dir
    V_diff <- var_p_haz_plogisN(mu, sigma, p_smaller) - var_dir
    E_diff^2 + V_diff^2
  }

  gamma_mu <- rep(0, D)
  sigma_gamma <- rep(0, D)
  for (i in 1:D) {
    theta_start <- if (i == 1) {
      c(-3, log(1))
    } else {
      c(gamma_mu[i - 1], log(sigma_gamma[i - 1]))
    }
    optim_res <- optim(theta_start, function(theta) {
      sum_sq_diff(theta, p_smaller = (i - 1) * e_dir)
    })
    gamma_mu[i] <- ifelse(optim_res$convergence == 0, optim_res$par[1], NA)
    sigma_gamma[i] <- ifelse(optim_res$convergence == 0, exp(optim_res$par[2]), NA)
  }
  data.frame(gamma_mu = gamma_mu, sigma_gamma = sigma_gamma)
}

## -----------------------
## iid Gamma prior
## -----------------------

# This is an alternative prior from a model by Höhle & an der Heiden
# Code was adjusted to obviate the requirement of the surveillance package

get_prior_gamma_iid <- function(prep_data_complete, now, start_date) {
  exp_var <- prep_data_complete %>%
    filter(event1_date >= start_date, event1_date <= now) %>%
    group_by(event1_date) %>%
    count() %>%
    ungroup() %>%
    complete(event1_date = seq.Date(start_date, now, by = "1 day"), fill = list(n = 0)) %>%
    summarize(mean = mean(n), var = var(n), .groups = "drop")

  exp <- exp_var$mean
  var <- exp_var$var

  dslnex <- function(x, aim = c(exp, var)) {
    y <- numeric(2)
    y[1] <- (x[1] * x[2]) - aim[1]
    y[2] <- (x[1] * x[2]^2) - aim[2]
    y
  }

  x_start <- rep(sqrt(exp), 2)
  res <- nleqslv(x_start, dslnex, control = list(btol = .01))
  beta.lambda <- 1 / res$x[2]
  alpha.lambda <- res$x[1]

  return(data.frame(alpha.lambda = alpha.lambda, beta.lambda = beta.lambda))
}
