# Utils for model fitting ----
## Sampling ----
define_sampling <- function(iter_warmup = 1000,
                            iter_sampling = 1000,
                            adapt_delta = 0.8,
                            max_treedepth = 10,
                            step_size = NULL,
                            chains = 4,
                            parallel_chains = NULL,
                            threads_per_chain = 1,
                            n_shards = NULL,
                            inits_function = get_inits,
                            seed = 42,
                            refresh = 200,
                            show_messages = T,
                            temp_save = T) {
  sampling_def <- list()
  sampling_def[["iter_warmup"]] <- iter_warmup
  sampling_def[["iter_sampling"]] <- iter_sampling
  sampling_def[["adapt_delta"]] <- adapt_delta
  sampling_def[["step_size"]] <- step_size
  sampling_def[["max_treedepth"]] <- max_treedepth
  sampling_def[["chains"]] <- chains
  sampling_def[["parallel_chains"]] <- ifelse(is.null(parallel_chains), chains, parallel_chains)
  sampling_def[["threads_per_chain"]] <- threads_per_chain
  sampling_def[["n_shards"]] <- ifelse(is.null(n_shards), threads_per_chain, n_shards)
  sampling_def[["inits_function"]] <- inits_function
  sampling_def[["seed"]] <- seed
  sampling_def[["refresh"]] <- refresh
  sampling_def[["show_messages"]] <- show_messages
  sampling_def[["temp_save"]] <- temp_save
  return(sampling_def)
}

define_output <- function(fit = TRUE,
                          model_type = TRUE,
                          job_id = TRUE,
                          date_index = TRUE,
                          dataset_index = TRUE,
                          stan_data_list = TRUE,
                          stan_prep_info = TRUE,
                          priors = TRUE,
                          summary = TRUE,
                          diagnostics = TRUE,
                          diagnostic_summary = TRUE,
                          profiles = TRUE,
                          output = TRUE,
                          delays = FALSE,
                          posterior_nowcast = FALSE,
                          posterior_R = FALSE,
                          overwrite = TRUE) {
  output_def <- c()

  args <- as.list(environment())
  for (a in names(args)) {
    if (a != "output_def") {
      if (args[[a]]) output_def <- c(output_def, a)
    }
  }
  return(output_def)
}

## Inits ----
random_inits <- function(stan_data_list, model_type) {
  return(NULL)
}

default_inits <- function(stan_data_list, model_type) {
  known_by_reference_date <- rowSums(stan_data_list$reported_known)
  if (sum(stan_data_list$reported_unknown) > 0) {
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
    alpha_emp <- rep(1 - 1e-4, 7)
  }
  alpha_emp_logit <- qlogis(alpha_emp)

  # get crude estimate of baseline hazard, assume beta and eta are zero
  gamma_emp <- rowMeans(sapply(1:7, function(i) stan_data_list$reported_known[i, ] / sum(stan_data_list$reported_known[i, ])))
  gamma_emp_haz <- get_hazard_from_p(gamma_emp)
  gamma_emp_haz_reduced <- rep(0, stan_data_list$n_delays)
  for (i in 1:stan_data_list$n_delays) {
    gamma_emp_haz_reduced[i] <- mean(gamma_emp_haz[stan_data_list$delay_idx == i])
  }

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
  last_unbiased <- length(I_by_reference_date) - stan_data_list$D
  I_by_reference_date <- c(I_by_reference_date[1:last_unbiased], rep(I_by_reference_date[last_unbiased], stan_data_list$D))
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

  init_function <- function() {
    inits <- list()

    # use an empirical estimate of alpha and assume it remains roughly unchanged over time
    inits[["alpha_logit_start"]] <- max(alpha_emp_logit, -10)
    inits[["alpha_logit_noise"]] <- rep(0, stan_data_list$T)
    inits[["alpha_logit_sd"]] <- max(stan_data_list$alpha_logit_sd_prior_mu + stan_data_list$alpha_logit_sd_prior_sd / 2, 1e-4)

    inits[["gamma"]] <- gamma_emp_haz_reduced + 1e-4
    inits[["beta"]] <- rep(1e-4, stan_data_list$n_beta)
    inits[["beta_sd"]] <- 0.5
    inits[["eta"]] <- rep(1e-4, stan_data_list$n_eta)

    inits[["xi_negbinom"]] <- max(stan_data_list$xi_negbinom_prior_mu + stan_data_list$xi_negbinom_prior_sd / 2, 1e-4)

    if (model_type %in% c("base", "base_old", "nowcast_imputed")) {
      inits[["lambda_log_start_values"]] <- c(log(stan_data_list$expected_cases_start), rep(0, 1 + stan_data_list$ets_diff))
      inits[["lambda_log_noise"]] <- rep(0, stan_data_list$T + stan_data_list$n_lambda_pre - 1 - stan_data_list$ets_diff) # rnorm(n = stan_data_list$T + stan_data_list$n_lambda_pre - 1)
      inits[["lambda_log_sd"]] <- max(stan_data_list$lambda_log_sd_prior_mu + stan_data_list$lambda_log_sd_prior_sd / 2, 1e-4) # rtnorm(n = 1, stan_data_list$lambda_log_sd_prior_mu, sd = stan_data_list$lambda_log_sd_prior_sd / 4, a = 0) # truncated
    }

    inits[["ets_alpha"]] <- 1 - 1e-4 # rbeta(n = 1, stan_data_list$ets_alpha_prior_alpha, stan_data_list$ets_alpha_prior_beta)
    inits[["ets_beta"]] <- 1 - 1e-4 # rbeta(n = 1, stan_data_list$ets_beta_prior_alpha, stan_data_list$ets_beta_prior_beta)
    inits[["ets_phi"]] <- 1 - 1e-4 # rbeta(n = 1, stan_data_list$ets_phi_prior_alpha, stan_data_list$ets_phi_prior_beta)

    inits[["R_level_start"]] <- R_emp[stan_data_list$max_gen + 1] # rtnorm(n = 1, mean = stan_data_list$R_level_start_prior_mu, sd = stan_data_list$R_level_start_prior_sd / 4, a = 0) # truncated
    inits[["R_trend_start"]] <- 1e-4 # rnorm(n = 1, mean = stan_data_list$R_trend_start_prior_mu, sd = stan_data_list$R_trend_start_prior_sd / 4)
    inits[["R_sd"]] <- 1
    inits[["R_noise"]] <- rep(0, stan_data_list$L + stan_data_list$n_lambda_pre + stan_data_list$T - stan_data_list$max_gen - 1) # c(0, diff(R_emp))[(stan_data_list$max_gen + 2):(stan_data_list$L + stan_data_list$n_lambda_pre + stan_data_list$T)]
    # use initial expected cases as proxy for infections / expected infections
    inits[["iota_log_ar_start"]] <- log(I_by_reference_date[1])
    inits[["iota_log_ar_sd"]] <- 1
    inits[["iota_log_ar_noise"]] <- c(0, diff(log(I_by_reference_date))[1:(stan_data_list$max_gen - 1)])
    inits[["I"]] <- I_by_reference_date + 0.1
    return(inits)
  }

  return(init_function)
}

## Stan utils ----
update_compiled_stanmodel <- function(model_def, force_recompile = FALSE) {
  n_models <- length(model_def$model_name)

  model_path_list <- lapply(1:n_models, function(i) here::here("code", model_def$model_folder, model_def$model_name[[i]]))
  include_paths_list <- lapply(1:n_models, function(i) here::here("code", model_def$model_folder)) # identical

  if (!model_def$profile) {
    stan_no_profiling_list <- lapply(1:n_models, function(i) write_stan_files_no_profiling(model_path_list[[i]], include_paths_list[[i]]))
    model_path_list <- lapply(1:n_models, function(i) stan_no_profiling_list[[i]]$model)
    include_paths_list <- lapply(1:n_models, function(i) stan_no_profiling_list[[i]]$include_paths)
  }

  stanmodel_list <- lapply(1:n_models, function(i) {
    cmdstan_model(model_path_list[[i]],
      include_paths = include_paths_list[[i]],
      cpp_options = list(stan_threads = model_def$threads),
      force_recompile = force_recompile
    )
  })

  model_def[["get_stan_model"]] <- lapply(1:n_models, function(i) {
    function() {
      return(stanmodel_list[[i]])
    }
  })

  return(model_def)
}
