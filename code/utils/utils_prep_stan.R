## ----------------------------------------------------------------
##                  Stan data list preparation                   -
## ----------------------------------------------------------------

prepare_data_list <- function(prep_data_complete,
                              prep_data_missing,
                              now,
                              start_date,
                              holidays = NULL,
                              D = 21,
                              delay_resolution = NULL,
                              model_type = "base",
                              gamma_prior_input = NULL,
                              gamma_prior_dirichlet_kappa = 1) {

  # --------------------------------------------
  # Reporting triangle
  reporting_matrices <- prepare_reporting_data(prep_data_complete, prep_data_missing, now, start_date, D)

  # --------------------------------------------
  # Delay resolution
  if (is.null(delay_resolution)) {
    delay_idx <- seq(1, D + 2)
  } else {
    stopifnot(delay_resolution[length(delay_resolution)] == Inf) # last delay resolution must be Inf
    delay_resolution <- delay_resolution[-length(delay_resolution)] # remove last delay resolution (we do not model it)
    stopifnot(sum(delay_resolution) == D) # ensure consistency with D
    delay_idx <- rep(seq(1, length(delay_resolution)), times = delay_resolution)
    delay_idx <- c(delay_idx, max(delay_idx) + 1, max(delay_idx) + 2) # add final elements for D and beyond
  }

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
  get_wdays <- function(nc_dates_t, D){
    wdays <- lubridate::wday(nc_dates_t + 0:(D - 1), label = TRUE)
    if(!is.null(holidays)){
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
  if (!is.null(delay_resolution)) {
    W <- reduce_design_matrix(W, delay_idx)
  }

  # ---------------------------
  # Prior for gamma
  if (!is.null(gamma_prior_input)) {
    stopifnot(nrow(gamma_prior_input) == D)
    gamma_prior_kappa <- get_prior_gamma(gd.prior.kappa = gamma_prior_dirichlet_kappa, D = D)
  } else {
    gamma_prior_kappa <- gamma_prior_input
  }

  gamma_prior_kappa_df <- data.frame(
    gamma_idx = delay_idx[1:(length(delay_idx) - 2)],
    gamma_mu = gamma_prior_kappa$gamma_mu,
    sigma_gamma = gamma_prior_kappa$sigma_gamma
  ) %>%
    group_by(gamma_idx) %>%
    summarize(across(everything(), mean), .groups = "drop")

  # ---------------------------
  # Mean of prior for initial iota
  iota_initial_mean <- rep(30, maxGen)
  iota_initial_sd <- rep(10, maxGen)

  #-------------------------------------------
  # Combine all data/preprocessed objects into list
  stan_data_list <- list(
    T = T + 1,
    D = D,
    n_delays = max(delay_idx) - 2,
    delay_idx = delay_idx,
    n_beta = length(ddChangepoint),
    n_eta = length(wdays),
    reported_known = reporting_matrices[["reported_known"]],
    reported_unknown = reporting_matrices[["reported_unknown"]],
    Z = Z,
    W = W,
    beta_mu = rep(0, length(ddChangepoint)),
    beta_sd = c(rep(0.01, length(ddChangepoint))),
    eta_mu = rep(0, length(wdays)),
    eta_sd = rep(0.5, length(wdays)),
    gamma_mu = gamma_prior_kappa_df$gamma_mu,
    gamma_sd = gamma_prior_kappa_df$sigma_gamma
  )

  if (model_type == "latent" | model_type == "renewal") {
    stan_data_list[["L"]] <- L
    stan_data_list[["latent_delay_dist"]] <- latent_delay_dist
  }

  if (model_type == "renewal") {
    stan_data_list[["max_gen"]] <- maxGen
    stan_data_list[["generation_time_dist"]] <- generation_time_dist
    stan_data_list[["iota_initial_mean"]] <- iota_initial_mean
    stan_data_list[["iota_initial_sd"]] <- iota_initial_sd
  }

  additional_info <- list(
    all_dates = all_dates,
    gamma_prior_kappa = gamma_prior_kappa
  )

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
