define_priors <- function(model_def,
                          dirichlet_prior = FALSE,
                          dirichlet_kappa = 1,
                          gamma_prior_precomputed_dir = NULL,
                          reporting_proportion = 1,
                          ...) {
  prior_def <- list()

  prior_def[["reporting_proportion"]] <- reporting_proportion

  if (dirichlet_prior) {
    prior_def[["dirichlet_kappa"]] <- dirichlet_kappa
    if (is.null(gamma_prior_precomputed_dir)) {
      gamma_prior_precomputed_dir <- here::here("code", model_def$model_folder, "priors_precomputed")
    }
    # Check if precomputed gamma (baseline reporting hazard) prior exists
    if (!dir.exists(gamma_prior_precomputed_dir)) dir.create(gamma_prior_precomputed_dir)
    gamma_prior_filepath <- here::here(gamma_prior_precomputed_dir, paste0("gamma_prior_k", dirichlet_kappa, "_D", model_def$D, ".rds"))
    prior_def[["gamma_prior_precomputed_path"]] <- gamma_prior_filepath
    if (file.exists(gamma_prior_filepath)) {
      gamma_prior_kappa <- readRDS(gamma_prior_filepath)
      stopifnot(nrow(gamma_prior_kappa) == model_def$D)
    } else {
      # compute prior
      gamma_prior_kappa <- get_prior_gamma(gd.prior.kappa = dirichlet_kappa, D = model_def[["D"]])
    }
    prior_def[["gamma_prior_kappa"]] <- gamma_prior_kappa

    gamma_prior_kappa_df <- data.frame(
      gamma_idx = model_def$delay_idx[1:(length(model_def$delay_idx) - 2)],
      gamma_mu = gamma_prior_kappa$gamma_mu,
      sigma_gamma = gamma_prior_kappa$sigma_gamma
    ) %>%
      group_by(gamma_idx) %>%
      summarize(across(everything(), mean), .groups = "drop")
  } else {
    # piecewise exponential prior for delay with mean such that P(Delay=maxDelay)=0.01
    mean_haz_logit <- qlogis(1 - (0.01)^(1 / model_def$D))
    # 95% upper bound for hazard such that P(Delay<=1)=0.90
    max_haz_logit <- qlogis(0.90)
    haz_sd_logit <- (max_haz_logit - mean_haz_logit) / 2

    gamma_prior_kappa_df <- data.frame(
      gamma_mu = rep(mean_haz_logit, model_def$D),
      sigma_gamma = rep(haz_sd_logit, model_def$D)
    )
  }

  additional_priors <- list(...)

  prior_def[["get_priors"]] <- function(stan_data_list) {
    priors <- list()

    # assumptions behind priors in comments
    priors[["alpha_logit_sd"]] <- get_prior("alpha_logit_sd", mu = 0, sd = 0.5) # odds ratio for missingness changes at max by a factor of 2.7 within one day
    priors[["alpha_logit_start"]] <- get_prior("alpha_logit_start", mu = 0, sd = 2) # alpha is between 2% and 98%
    priors[["beta"]] <- get_prior("beta", mu = rep(0, stan_data_list$n_beta), sd = c(rep(0.1, stan_data_list$n_beta))) # odds ratio for reporting hazard changes at max by a factor of 4 within one week
    if (model_def$delay_changepoint == "rw") {
      priors[["beta_random"]] <- TRUE
      priors[["beta_sd"]] <- get_prior("beta_sd", mu = 0.35, sd = 0.175) # odds ratio for reporting hazard changes at max by a factor of 4 within one week
    } else {
      priors[["beta_random"]] <- FALSE
      priors[["beta_sd"]] <- get_prior("beta_sd", mu = numeric(0), sd = numeric(0))
    }
    priors[["eta"]] <- get_prior("eta", mu = rep(0, stan_data_list$n_eta), sd = c(rep(0.5, stan_data_list$n_eta))) # odds ratio for reporting hazard changes at max by a factor of 2.7 for one weekday
    priors[["gamma"]] <- get_prior("gamma", mu = gamma_prior_kappa_df$gamma_mu, sd = gamma_prior_kappa_df$sigma_gamma) # see above

    priors[["xi_negbinom"]] <- get_prior("xi_negbinom", mu = 0, sd = 1) # generic prior on 1/sqrt(phi), stan recommendations

    priors[["ets_alpha_fixed"]] <- -1
    # Mean = 0.5, hence the last 7 days (including today) cover over 99% of weight
    # --> Lower 95% interval = 0.25, hence the last 20 days (including today) cover over 99% of weight
    # --> Upper 95% interval = 0.75, hence the last 3 days (including today) cover over 99% of weight
    # ets_total_weight <- function(alpha, days) return(alpha * sum((1-alpha)^(0:(days-1))))}
    priors[["ets_alpha"]] <- get_prior("ets_alpha", alpha = 5, beta = 5)
    priors[["ets_beta_fixed"]] <- -1
    priors[["ets_beta"]] <- get_prior("ets_beta", alpha = 5, beta = 5)
    priors[["ets_phi_fixed"]] <- -1
    priors[["ets_phi"]] <- get_prior("ets_phi", alpha = 50, beta = 5) # Mean = 0.9, 95% interval 0.82-0.97, see recommended values for ETS
    # if mean = 0.9, then trend only half as strong after 7 days (0.9^7)

    if (model_def$model_type %in% c("base", "base_old", "nowcast_imputed", "impute_and_nowcast", "impute_parametric_and_nowcast")) {
      priors[["lambda_log_sd"]] <- get_prior("lambda_log_sd", mu = 0, sd = 0.05)
      priors[["lambda_log_level_start"]] <- get_prior("lambda_log_level_start", mu = log(stan_data_list$expected_cases_start + 0.1), sd = 0.5)
      priors[["lambda_log_trend_start"]] <- get_prior("lambda_log_trend_start", mu = 0, sd = 0.05)
      priors[["lambda_log_2nd_trend_start"]] <- get_prior("lambda_log_2nd_trend_start", mu = 0, sd = 0.02)
    }

    priors[["R_sd"]] <- get_prior("R_sd", mu = 0, sd = 0.1) # note that this is (approximately) on the absolute scale: expect R to change at max by 1.4 within one week
    priors[["R_level_start"]] <- get_prior("R_level_start", mu = 1, sd = 0.8) # default R = 1, 95% interval between 2.6 and close to zero
    priors[["R_trend_start"]] <- get_prior("R_trend_start", mu = 0, sd = 0.1) # mean = no trend , 95% interval is 1.4 per day
    priors[["iota_log_ar_start"]] <- get_prior("iota_log_ar_start", mu = log(stan_data_list$expected_cases_start / reporting_proportion), sd = 0.2) # slightly informed by data (expected cases), 95% interval: *0.67 -- *1.5
    priors[["iota_log_ar_sd"]] <- get_prior("iota_log_ar_sd", mu = 0.05, sd = 0.025) # log daily growth rate of up to 0.2, i.e. roughly relative changes of 20% (same for negative growth rates)

    # add further priors, potentially overwriting existing
    if (length(additional_priors) > 0) {
      for (i in 1:length(additional_priors)) {
        priors[[names(additional_priors)[i]]] <- additional_priors[[i]]
      }
    }
    for (ets_param in c("ets_alpha", "ets_beta", "ets_phi")) {
      if (paste0(ets_param, "_fixed") %in% names(priors)) {
        if (priors[[paste0(ets_param, "_fixed")]] >= 0) {
          priors[[ets_param]] <- get_prior(ets_param, alpha = numeric(0), beta = numeric(0))
        }
      }
    }
    if (stan_data_list$ets_diff == 0) {
      priors[["lambda_log_2nd_trend_start"]] <- get_prior("lambda_log_2nd_trend_start", mu = numeric(0), sd = numeric(0))
    }

    return(priors)
  }
  return(prior_def)
}

##                Prior for the delay distribution

##  Dirichlet gamma prior

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

## iid Gamma prior
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
