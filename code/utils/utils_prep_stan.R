## ----------------------------------------------------------------
##                  Stan data list preparation                   -
## ----------------------------------------------------------------

get_inits <- function(stan_data_list, model_type) {
  return(
    function() {
      inits <- list()

      inits[["alpha_logit_start"]] <- rnorm(n = 1, mean = 0, sd = 2)
      inits[["alpha_logit_raw"]] <- rep(rnorm(n = 1), stan_data_list[[T]])
      inits[["alpha_logit_sd"]] <- rhnorm(n = 1, sigma = 0.5) # half-normal

      inits[["gamma"]] <- rnorm(n = stan_data_list$n_delays, mean = stan_data_list$gamma_prior_mu, sd = stan_data_list$gamma_prior_sd)
      inits[["beta"]] <- rnorm(n = stan_data_list$n_beta, mean = stan_data_list$beta_prior_mu, sd = stan_data_list$beta_prior_sd)
      inits[["eta"]] <- rnorm(n = stan_data_list$n_eta, mean = stan_data_list$eta_prior_mu, sd = stan_data_list$eta_prior_sd)

      inits[["xi_negbinom"]] <- rnorm(n = 1, mean = stan_data_list$xi_negbinom_prior_mu, sd = stan_data_list$xi_negbinom_prior_sd)
      
      if (model_type == "base") {
        inits[["lambda_log_start"]] <- rnorm(n = 1, mean = log(stan_data_list$expected_cases_start), sd = 0.2)
        inits[["lambda_log_raw"]] <- rnorm(n = stan_data_list[[T]])
        inits[["lambda_log_sd"]] <- rhnorm(n = 1, sigma = 0.5) # half-normal
      } else if (model_type == "latent") {
        # use initial expected cases as proxy for initial expected infections
        inits[["iota_log_start"]] <- rnorm(n = 1, mean = log(stan_data_list$expected_cases_start), sd = 0.4)
        inits[["iota_log_raw"]] <- rnorm(n = stan_data_list[[T]])
        inits[["iota_log_sd"]] <- rhnorm(n = 1, sigma = 0.5) # half-normal
      } else if (model_type == "renewal") {
        inits[["R_ets_alpha"]] <- rbeta(n = 1, 4, 4)
        inits[["R_ets_beta"]] <- rbeta(n = 1, 4, 4)
        inits[["R_ets_phi"]] <- rbeta(n = 1, 50, 5)
        inits[["R_level_start"]] <- rnorm(n = 1, mean = log(1), 0.8)
        inits[["R_trend_start"]] <- rnorm(n = 1, mean = 0, 0.2)
        inits[["R_raw"]] <- rnorm(n = stan_data_list[[T]])
        inits[["R_sd"]] <- rhnorm(n = 1, sigma = 0.3) # half-normal

        # use initial expected cases as proxy for infections / expected infections
        inits[["iota_initial"]] <- rtnorm(n = stan_data_list$max_gen, mean = stan_data_list$expected_cases_start, sd = 0.5, a = 0)
        I_length <- stan_data_list$max_gen + stan_data_list$L + stan_data_list$D + stan_data_list$T
        inits[["I"]] <- rtnorm(n = I_length, mean = stan_data_list$expected_cases_start, sd = 0.6, a = 0)
      }

      return(inits)
    }
  )
}

get_prior <- function(variable, ...){
  prior_list <- list()
  prior_values <- list(...)
  for(v in names(prior_values)){
    prior_list[[paste0(variable,"_prior_",v)]] <- prior_values[[v]]
  }
  return(prior_list)
}

define_priors <- function(model_def,
                          dirichlet_kappa = 1,
                          gamma_prior_precomputed_dir = here::here("code", "models", "priors_precomputed")) {
  prior_def <- list()
  prior_def[["dirichlet_kappa"]] <- dirichlet_kappa
  
  # Check if precomputed gamma (baseline reporting hazard) prior exists
  if (!dir.exists(gamma_prior_precomputed_dir)) dir.create(gamma_prior_precomputed_dir)
  gamma_prior_filepath <- here::here(gamma_prior_precomputed_dir, paste0("gamma_prior_k", dirichlet_kappa, "_D", model_def$D, ".rds"))
  prior_def[["gamma_prior_precomputed_path"]] <- gamma_prior_filepath
  if (file.exists(gamma_prior_filepath)) {
    gamma_prior_kappa <- readRDS(gamma_prior_filepath)
    stopifnot(nrow(gamma_prior_kappa) == model_def$D)
  } else {
    # precompute prior
    gamma_prior_kappa <- get_prior_gamma(gd.prior.kappa = dirichlet_kappa, D = D)
  }
  prior_def[["gamma_prior_kappa"]] <- gamma_prior_kappa
  
  gamma_prior_kappa_df <- data.frame(
    gamma_idx = model_def$delay_idx[1:(length(model_def$delay_idx) - 2)],
    gamma_mu = gamma_prior_kappa$gamma_mu,
    sigma_gamma = gamma_prior_kappa$sigma_gamma
  ) %>%
    group_by(gamma_idx) %>%
    summarize(across(everything(), mean), .groups = "drop")
  
  prior_def[["get_priors"]] <- function(stan_data_list){
    priors <- list()
    
    priors[["alpha_logit_sd"]] = get_prior("alpha_logit_sd", mu = 0, sd = 0.5)
    priors[["alpha_logit_start"]] = get_prior("alpha_logit_start", mu = 0, sd = 2)
    priors[["beta"]] = get_prior("beta", mu = rep(0, stan_data_list$n_beta), sd = c(rep(0.01, stan_data_list$n_beta)))
    priors[["eta"]] = get_prior("eta", mu = rep(0, stan_data_list$n_eta), sd = c(rep(0.5, stan_data_list$n_eta)))
    priors[["gamma"]] = get_prior("gamma", mu = gamma_prior_kappa_df$gamma_mu, sd = gamma_prior_kappa_df$sigma_gamma)
    
    priors[["xi_negbinom"]] = get_prior("xi_negbinom", mu = 0, sd = 1)
    
    if (model_def$model_type == "base") {
      priors[["lambda_log_sd"]] <- get_prior("lambda_log_sd", mu = 0, sd = 0.5)
      priors[["lambda_log_start"]] <- get_prior("lambda_log_start", mu = 0, sd = 12)
    }
    
    if (model_def$model_type == "latent") {
      priors[["iota_log_sd"]] <- get_prior("iota_log_sd", mu = 0, sd = 0.5)
      priors[["iota_log_start"]] <- get_prior("iota_log_start", mu = 0, sd = 12)
    }
    
    if (model_def$model_type == "renewal") {
      priors[["R_ets_alpha"]] <- get_prior("R_ets_alpha", alpha = 70, beta = 30) # the last three days have >5% impact
      priors[["R_ets_beta"]] <- get_prior("R_ets_beta", alpha = 30, beta = 70) # the last six days have >5% impact
      priors[["R_ets_phi"]] <- get_prior("R_ets_phi", alpha = 50, beta = 5)
      priors[["R_sd"]] <- get_prior("R_sd", mu = 0, sd = 0.3)
      priors[["R_level_start"]] <- get_prior("R_level_start", mu = log(1), sd = log(5))
      priors[["R_trend_start"]] <- get_prior("R_trend_start", mu = log(1), sd = log(1.1))
      priors[["iota_initial"]] <- get_prior("iota_initial", mu = rep(30, stan_data_list$max_gen), sd = rep(10, stan_data_list$max_gen))
    }
    
    return(priors)
  }
  
  return(prior_def)
  
}

prepare_data_list <- function(prep_data_complete,
                              prep_data_missing,
                              now,
                              start_date,
                              holidays = NULL,
                              D,
                              delay_idx,
                              model_type = "base",
                              get_priors) {

  # --------------------------------------------
  # Reporting triangle
  reporting_matrices <- prepare_reporting_data(prep_data_complete, prep_data_missing, now, start_date, D)

  # estimate initial expected cases based on reported known + unknown assuming a flat delay distribution
  expected_cases_start <- sum(reporting_matrices[["reported_known"]][1, ]) +
    mean(reporting_matrices[["reported_unknown"]][1:D])

  # ---------------------------
  # Covariates for reporting delay
  all_dates <- seq(start_date, now, by = "1 day")
  tswp <- seq(start_date - D, now, by = "1 day") # time series including the days of the modeling phase
  T <- length(all_dates) - 1

  ## 1) Segmented linear model for change over time
  ddChangepoint <- rev(seq(now, start_date, by = "-2 weeks"))
  Z_idx <- rev(rep(length(ddChangepoint):1, each = 14)[1:(1 + T)])
  Z <- array(NA, dim = c(T + 1 + D, length(ddChangepoint)), dimnames = list(as.character(tswp), 1:length(ddChangepoint)))
  for (t_pre in 1:D) {
    Z[t_pre, ] <- rep(0, length(ddChangepoint)) # we don't have a changepoint model before the inference phase, but fill the space up with 0s to have the same dimensionality as the other matrices
  }
  for (t in 0:T) {
    covariates_ind_days <- c(0, as.numeric(all_dates <= (all_dates[t + 1]))[-1])
    Z[D + 1 + t, ] <- tapply(covariates_ind_days, Z_idx, sum) # add up counts in each changepoint interval
  }
  # scale columns by standard deviation for better efficiency
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
  if (!all(delay_idx==seq(1, D + 2))) {
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
  
  # observed data
  stan_data_list[["reported_known"]] <- reporting_matrices[["reported_known"]]
  stan_data_list[["reported_unknown"]] <- reporting_matrices[["reported_unknown"]]
  stan_data_list[["expected_cases_start"]] <- expected_cases_start
  
  # reporting model
  stan_data_list[["n_delays"]] <- max(delay_idx) - 2
  stan_data_list[["delay_idx"]] <- delay_idx
  stan_data_list[["n_beta"]] <- length(ddChangepoint)
  stan_data_list[["n_eta"]] <- length(wdays)

  # latent process details
  if (model_type == "latent" | model_type == "renewal") {
    stan_data_list[["L"]] <- L
    stan_data_list[["latent_delay_dist"]] <- latent_delay_dist
  }

  # renewal process details
  if (model_type == "renewal") {
    stan_data_list[["max_gen"]] <- maxGen
    stan_data_list[["generation_time_dist"]] <- generation_time_dist
  }
  
  # all dates
  additional_info[["all_dates"]] <- all_dates
  
  # priors
  priors <- get_priors(stan_data_list)
  additional_info[["priors"]] <- priors
  stan_data_list <- append(stan_data_list, flatten(priors))
  
  # inits
  additional_info[["inits"]] <- get_inits(stan_data_list, model_type)

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
    summarize(n = sum(n), .groups = "drop") %>%
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
  unobserved_matrix <- prep_data_missing %>%
    group_by(event2_date) %>%
    summarize(n = sum(n), .groups = "keep") %>%
    filter(event2_date >= start_date, event2_date <= now) %>%
    ungroup() %>%
    complete(event2_date = seq.Date(start_date, now, by = "1 day"), fill = list(n = 0)) %>%
    pull(n)

  return(list(reported_known = reporting_matrix, reported_unknown = unobserved_matrix))
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
