#' Helper function to define a prior for stan
#'
#' @param variable Name of the variable for which the prior is defined
#' @param ... Named arguments representing the arguments of the distribution
#'   function
#'
#' @return
get_prior <- function(variable, ...) {
  prior_list <- list()
  prior_values <- list(...)
  for (v in names(prior_values)) {
    prior_list[[paste0(variable, "_prior_", v)]] <- prior_values[[v]]
  }
  return(prior_list)
}

#' Get (default) priors for nowcasting models
#'
#' @param model_def Model information for context, see `define_model()`
#' @param reporting_proportion Proportion of infections which become reported
#'   cases. 1 means full reporting (100% of infections are reported).
#' @param latent_delay_dist A `vector` representing a discrete distribution of
#'   the delay between latent (e.g. infections) and observed (e.g. symptom
#'   onsets) events, e.g. the incubation period.
#' @param generation_time_dist A `vector` representing a discrete generation
#'   time distribution.
#' @param dirichlet Should an (approximate) Dirichlet prior be used for the
#'   baseline reporting delay?
#' @param ... Further user-defined priors, can also be used to overwrite
#'   defaults
#'
#' @return A `list` with prior information, in particular an element
#'   `get_priors`, which is a function that takes a stan_data_list as input and
#'   returns suitable priors for all parameters.
define_priors <- function(model_def,
                          reporting_proportion = 1,
                          latent_delay_dist = NULL,
                          generation_time_dist = NULL,
                          dirichlet = FALSE,
                          ...) {
  prior_def <- list()

  prior_def[["reporting_proportion"]] <- reporting_proportion
  prior_def[["latent_delay_dist"]] <- latent_delay_dist
  prior_def[["generation_time_dist"]] <- generation_time_dist
  
  gamma_prior_df <- get_prior_gamma(model_def, dirichlet = dirichlet)

  additional_priors <- list(...)

  prior_def[["get_priors"]] <- function(stan_data_list) {
    priors <- list()

    # Intercept for probability of known onset date
    # alpha is between 2% and 98%
    priors[["alpha_logit_start"]] <- get_prior("alpha_logit_start", mu = 0, sd = 2)

    # Standard deviation of random walk on $\text{logit}
    # odds for missingness changes at max by a factor of 1/3 or 3 within one day
    priors[["alpha_logit_sd"]] <- get_prior("alpha_logit_sd", mu = 0, sd = 0.5)

    # Linear trend of the logit reporting hazard during week $i$
    # odds for reporting hazard changes at max by a factor of 4 within one week
    priors[["beta"]] <- get_prior(
      "beta",
      mu = rep(0, stan_data_list$n_beta),
      sd = c(rep(0.1, stan_data_list$n_beta))
    )
    if (model_def$delay_changepoint == "rw") {
      priors[["beta_random"]] <- TRUE
      # odds for reporting hazard changes at max by a factor of 4 within one week
      priors[["beta_sd"]] <- get_prior("beta_sd", mu = 0.35, sd = 0.175)
    } else {
      priors[["beta_random"]] <- FALSE
      priors[["beta_sd"]] <- get_prior("beta_sd", mu = numeric(0), sd = numeric(0))
    }

    # Effect of weekday $i$ on logit reporting hazard
    # odds for reporting hazard changes at max by a factor of 2.7 for one weekday
    priors[["eta"]] <- get_prior(
      "eta",
      mu = rep(0, stan_data_list$n_eta),
      sd = c(rep(0.75, stan_data_list$n_eta))
    )

    # Baseline logit reporting hazard for a delay of $d$ days
    priors[["gamma"]] <- get_prior(
      "gamma",
      mu = gamma_prior_df$gamma_mu,
      sd = gamma_prior_df$sigma_gamma
    )

    # Overdispersion parameter for case counts
    # generic prior on 1/sqrt(phi), stan dev team recommendations
    priors[["xi_negbinom"]] <- get_prior("xi_negbinom", mu = 0, sd = 1)

    # Exponential smoothing
    # Smoothing parameter for level
    # Mean = 0.5 => last 7 days (including today) cover over 99% of weight
    # Lower 95% interval = 0.25, => last 20 days (including today) cover over 99% of weight
    # Upper 95% interval = 0.75 => last 3 days (including today) cover over 99% of weight
    # ets_total_weight <- function(alpha, days) return(alpha * sum((1-alpha)^(0:(days-1))))}
    priors[["ets_alpha_fixed"]] <- -1
    priors[["ets_alpha"]] <- get_prior("ets_alpha", alpha = 5, beta = 5)

    # Smoothing parameter for trend
    # Same as for alpha
    priors[["ets_beta_fixed"]] <- -1
    priors[["ets_beta"]] <- get_prior("ets_beta", alpha = 5, beta = 5)

    # Dampening parameter for trend
    # Mean = 0.9, 95% interval 0.82-0.97, see recommended values for ETS
    # if 0.9, then trend only half as strong after 7 days (0.9^7)
    priors[["ets_phi_fixed"]] <- -1
    priors[["ets_phi"]] <- get_prior("ets_phi", alpha = 50, beta = 5)

    stepwise_models <- c(
      "impute_adjust",
      "adjust",
      "impute_then_adjust",
      "impute_independent_then_adjust",
      "impute_parametric_then_adjust"
    )

    if (model_def$model_type %in% stepwise_models) {
      # Standard deviation of random walk on $\log(\lambda_t)$
      # log daily growth rate of up to 0.2, i.e. roughly relative changes of
      # 20% (same for negative growth rates)
      priors[["lambda_log_sd"]] <- get_prior(
        "lambda_log_sd",
        mu = 0.05, sd = 0.025
      )

      # Intercept for $\log$ symptom onsets (and higher order diff intercepts)
      # slightly informed by data (expected cases) with x3 uncertainty margin
      priors[["lambda_log_level_start"]] <- get_prior(
        "lambda_log_level_start",
        mu = log(stan_data_list$expected_cases_start + 0.1),
        sd = 0.5
      )
      priors[["lambda_log_trend_start"]] <- get_prior(
        "lambda_log_trend_start",
        mu = 0,
        sd = 0.05
      )
      priors[["lambda_log_2nd_trend_start"]] <- get_prior(
        "lambda_log_2nd_trend_start",
        mu = 0,
        sd = 0.02
      )
    }

    # Intercept for $\log$ infections in seeding phase
    # slightly informed by data (expected infections) with x3 uncertainty margin
    priors[["iota_log_ar_start"]] <- get_prior(
      "iota_log_ar_start",
      mu = log(max(stan_data_list$expected_cases_start /
        reporting_proportion, 1e-4)),
      sd = 0.5
    )

    # Standard deviation of random walk on $\log(\iota_t)$
    # log daily growth rate of up to 0.2, i.e. relative changes of approx 20%
    # (same for negative growth rates)
    # Mote: if log sd is 0.1 at max, then log growth increment is 0.2 at max
    priors[["iota_log_ar_sd"]] <- get_prior("iota_log_ar_sd", mu = 0.05, sd = 0.025)

    # Intercept for effective reproduction number
    # default R = 1, 95% interval between 2.6 and close to zero
    priors[["R_level_start"]] <- get_prior("R_level_start", mu = 1, sd = 0.8)

    # Trend intercept for effective reproduction number
    # mean = no trend , 95% interval is 1.4 per day
    priors[["R_trend_start"]] <- get_prior("R_trend_start", mu = 0, sd = 0.1)

    # Standard deviation of random walk on $R_t$
    # note that this is (approximately) on the absolute scale:
    # we expect R to change at max by +-2.8 within one week
    priors[["R_sd"]] <- get_prior("R_sd", mu = 0, sd = 0.1)

    # Add user-defined priors, potentially overwriting existing
    if (length(additional_priors) > 0) {
      for (i in 1:length(additional_priors)) {
        priors[[names(additional_priors)[i]]] <- additional_priors[[i]]
      }
    }

    # Remove ETS priors in case of smoothing / dampening parameters are fixed
    for (ets_param in c("ets_alpha", "ets_beta", "ets_phi")) {
      if (paste0(ets_param, "_fixed") %in% names(priors)) {
        if (priors[[paste0(ets_param, "_fixed")]] >= 0) {
          priors[[ets_param]] <- get_prior(
            ets_param,
            alpha = numeric(0), beta = numeric(0)
          )
        }
      }
    }

    # Remove 2nd order trend if not required
    if (stan_data_list$ets_diff == 0) {
      priors[["lambda_log_2nd_trend_start"]] <- get_prior(
        "lambda_log_2nd_trend_start",
        mu = numeric(0), sd = numeric(0)
      )
    }
    return(priors)
  }
  return(prior_def)
}

#' Get a prior for the baseline reporting hazard
#'
#' @param model_def Model information for context, see `define_model()`
#' @param dirichlet Should an (approximated) Dirichlet prior be used? If FALSE,
#'   a piecewise exponential prior is used instead.
#' @param dirichlet_kappa Kappa parameter of the Dirichlet distribution
#' @param gamma_prior_precomputed_dir A directory where precomputed Dirichlet
#'   prior approximations are stored
#'
#' @return A `data.frame`with mu and sigma parameters for the normal
#'   distribution prior at each delay
get_prior_gamma <- function(model_def, dirichlet = FALSE, dirichlet_kappa = 1,
                            gamma_prior_precomputed_dir = NULL) {
  if (dirichlet) {
    if (is.null(gamma_prior_precomputed_dir)) {
      gamma_prior_precomputed_dir <- here::here(
        "code", "models", "priors_precomputed"
      )
    }

    # Check if precomputed gamma (baseline reporting hazard) prior exists
    if (!dir.exists(gamma_prior_precomputed_dir)) {
      dir.create(gamma_prior_precomputed_dir)
    }
    gamma_prior_filepath <- here::here(
      gamma_prior_precomputed_dir,
      paste0("gamma_prior_k", dirichlet_kappa, "_D", model_def$D, ".rds")
    )
    prior_def[["gamma_prior_precomputed_path"]] <- gamma_prior_filepath
    if (file.exists(gamma_prior_filepath)) {
      gamma_prior_kappa <- readRDS(gamma_prior_filepath)
      stopifnot(nrow(gamma_prior_kappa) == model_def$D)
    } else {
      # compute prior
      gamma_prior_kappa <- get_prior_gamma_dirichlet(
        gd.prior.kappa = dirichlet_kappa, D = model_def$D
      )
    }
    prior_def[["gamma_prior_kappa"]] <- gamma_prior_kappa

    gamma_prior_df <- data.frame(
      gamma_idx = model_def$delay_idx[1:(length(model_def$delay_idx) - 2)],
      gamma_mu = gamma_prior_kappa$gamma_mu,
      sigma_gamma = gamma_prior_kappa$sigma_gamma
    ) %>%
      group_by(gamma_idx) %>%
      summarize(across(everything(), mean), .groups = "drop")

    return(gamma_prior_df)
  } else {
    # piecewise exponential prior for delay
    # with mean such that P(Delay=maxDelay)=0.01
    mean_haz_logit <- qlogis(1 - (0.01)^(1 / model_def$D))

    # 95% upper bound for hazard such that P(Delay=0)=0.98
    max_haz_logit <- qlogis(0.98)
    haz_sd_logit <- (max_haz_logit - mean_haz_logit) / 2

    gamma_prior_df <- data.frame(
      gamma_idx = model_def$delay_idx[1:(length(model_def$delay_idx) - 2)],
      gamma_mu = rep(mean_haz_logit, model_def$D),
      sigma_gamma = rep(haz_sd_logit, model_def$D)
    ) %>%
      group_by(gamma_idx) %>%
      summarize(across(everything(), mean), .groups = "drop")

    return(gamma_prior_df)
  }
}

#' Approximate a Dirichlet prior for the baseline reporting hazard using
#' normally distributed priors
#'
#' The goal of the below functions is to find a good prior (through mu and sigma
#' of a normal distribution) for the gamma parameter of the model (baseline
#' logit-hazard for delay d). This is achieved by choosing the mu and sigma such
#' that the resulting delay distribution is as similar as possible to a Dirichlet
#' distribution where all days within the horizon D have the same probability.
#' The kappa parameter (related to the variance of the Dirichlet) can be chosen
#' by the user.
#'
#' Code taken from Günther et. al (2021), see
#' https://github.com/FelixGuenther/nc_covid19_bavaria, with minor improvement to
#' respect non-negativity of variance parameter
#'
#' @param gd.prior.kappa Kappa parameter of Dirichlet distribution
#' @param D Assumed maximum delay
#'
#' @return A `data.frame` with mu and sigma parameters for the normal
#'   distribution prior at each delay
get_prior_gamma_dirichlet <- function(gd.prior.kappa = 1, D = 20) {
  # Dirichlet gamma prior helper functions
  # Expectation of P(D>=d)*P(D=d|D>=d)
  expectation_p_haz_plogisN <- function(mu, sigma, p_smaller) {
    f <- function(x) {
      (1 - p_smaller) * plogis(x) * dnorm(x, mu, sigma)
    }
    int <- integrate(f, lower = -Inf, upper = Inf)
    int$value
  }

  # Variance of P(D>=d)*P(D=d|D>=d),
  # with P(D=d|D>=d)=1/(1+exp(-x)) and X~N(mu, sigma^2)
  var_p_haz_plogisN <- function(mu, sigma, p_smaller) {
    E <- expectation_p_haz_plogisN(mu, sigma, p_smaller)
    f <- function(x) {
      (((1 - p_smaller) * plogis(x)) - E)^2 * dnorm(x, mu, sigma)
    }
    int <- integrate(f, lower = -Inf, upper = Inf)
    int$value
  }

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

#' Get an iid Gamma prior for the baseline reporting hazard
#'
#' Alternative to the Dirichlet prior, taken from  Höhle & an der Heiden 2014,
#' code was adjusted to obviate the requirement of the surveillance package
get_prior_gamma_iid <- function(prep_data_complete, now, start_date) {
  exp_var <- prep_data_complete %>%
    filter(event1_date >= start_date, event1_date <= now) %>%
    group_by(event1_date) %>%
    count() %>%
    ungroup() %>%
    complete(
      event1_date = seq.Date(start_date, now, by = "1 day"),
      fill = list(n = 0)
    ) %>%
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
