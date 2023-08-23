#' Define parameters of the NUTS sampler in stan
#'
#' The arguments are mostly passed to (and thus identical to) `$sample()` method
#' of cmdstanr.
#'
#' @param iter_warmup (positive integer) The number of warmup iterations to run
#'   per chain. Note: in the CmdStan User's Guide this is referred to as
#'   num_warmup.
#' @param iter_sampling (positive integer) The number of post-warmup iterations
#'   to run per chain. Note: in the CmdStan User's Guide this is referred to as
#'   num_samples.
#' @param adapt_delta (real in (0,1)) The adaptation target acceptance
#'   statistic.
#' @param max_treedepth (positive integer) The maximum allowed tree depth for
#'   the NUTS engine. See the Tree Depth section of the CmdStan User's Guide for
#'   more details.
#' @param step_size (positive real) The initial step size for the discrete
#'   approximation to continuous Hamiltonian dynamics. This is further tuned
#'   during warmup.
#' @param chains (positive integer) The number of Markov chains to run. The
#'   default is 4.
#' @param parallel_chains (positive integer) The maximum number of MCMC chains
#'   to run in parallel. If parallel_chains is not specified then the default is
#'   to look for the option "mc.cores", which can be set for an entire R session
#'   by options(mc.cores=value). If the "mc.cores" option has not been set then
#'   the default is 1.
#' @param threads_per_chain (positive integer) If the model was compiled with
#'   threading support, the number of threads to use in parallelized sections
#'   within an MCMC chain (e.g., when using the Stan functions reduce_sum() or
#'   map_rect()). This is in contrast with parallel_chains, which specifies the
#'   number of chains to run in parallel. The actual number of CPU cores used is
#'   parallel_chains*threads_per_chain.
#' @param n_shards Number of shards to use when running within-chain
#'   parallelization. Currently not supported by most models.
#' @param inits_function A function that returns a single list with names
#'   corresponding to the parameters for which you are specifying initial
#'   values. The function can take no arguments or a single argument chain_id.
#'   For MCMC, if the function has argument chain_id it will be supplied with
#'   the chain id (from 1 to number of chains) when called to generate the
#'   initial values.
#' @param seed positive integer(s)) A seed for the (P)RNG to pass to CmdStan. In
#'   the case of multi-chain sampling the single seed will automatically be
#'   augmented by the the run (chain) ID so that each chain uses a different
#'   seed. The exception is the transformed data block, which defaults to using
#'   same seed for all chains so that the same data is generated for all chains
#'   if RNG functions are used. The only time seed should be specified as a
#'   vector (one element per chain) is if RNG functions are used in transformed
#'   data and the goal is to generate different data for each chain.
#' @param refresh (non-negative integer) The number of iterations between
#'   printed screen updates. If refresh = 0, only error messages will be
#'   printed.
#' @param show_messages (logical) When TRUE (the default), prints all
#'   informational messages, for example rejection of the current proposal.
#'   Disable if you wish to silence these messages, but this is not usually
#'   recommended unless you are very confident that the model is correct up to
#'   numerical error. If the messages are silenced then the $output() method of
#'   the resulting fit object can be used to display the silenced messages.
#' @param temp_save Should the fitted model be saved to a temporary file before
#'   postprocessing?
#'
#' @return A `list` with parameters for the NUTS sampler.
define_sampling <- function(iter_warmup = 1000,
                            iter_sampling = 1000,
                            adapt_delta = NULL,
                            max_treedepth = NULL,
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
  parallel_chains <- ifelse(is.null(parallel_chains), chains, parallel_chains)
  sampling_def[["parallel_chains"]] <- parallel_chains
  sampling_def[["threads_per_chain"]] <- threads_per_chain
  n_shards <- ifelse(is.null(n_shards), threads_per_chain, n_shards)
  sampling_def[["n_shards"]] <- n_shards
  sampling_def[["inits_function"]] <- inits_function
  sampling_def[["seed"]] <- seed
  sampling_def[["refresh"]] <- refresh
  sampling_def[["show_messages"]] <- show_messages
  sampling_def[["temp_save"]] <- temp_save
  if (sampling_def$threads_per_chain == 1) {
    sampling_def$threads_per_chain <- NULL
  }
  return(sampling_def)
}

#' Define which results of the model fitting should be returned
#'
#' @param fit The fitted model
#' @param model_type The type of model that was fitted
#' @param job_id Job identifier, if run on cluster
#' @param date_index Index of the nowcasting date at which the model was fitted
#'   (if part of a job array)
#' @param dataset_index Index of the dataset to which the model was fitted (if
#'   part of a job array)
#' @param stan_data_list List with all data past to stan
#' @param stan_prep_info Meta information from the preprocessing
#' @param priors The priors used by the model
#' @param summary Summary of nowcasting results from the model fit
#' @param diagnostics Detailed sampler diagnostics for each chain, as a
#'   `draws_array`.
#' @param diagnostic_summary Summaries of sampler diagnostics and warning
#'   messages.
#' @param profiles Profiling information about efficiency of stan code.
#' @param output The printed output from the chains.
#' @param delays Posterior distribution summary of delays for each combination
#'   of date of symptom onset and delay.
#' @param R_epiestim The reproduction number estimated using EpiEstim.
#' @param R_generative The reproduction number estimated using a generative
#'   renewal model.
#' @param posterior_nowcast Posterior samples of the nowcast for cases.
#' @param posterior_R Posterior samples of the reproduction number R.
#' @param overwrite If a model fit/result file with the same name already
#'   exists, should it be overwritten?
#'
#' @return A `vector` with all types of output/results that should be returned.
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
                          R_epiestim = TRUE,
                          R_generative = TRUE,
                          posterior_nowcast = FALSE,
                          posterior_R = FALSE,
                          overwrite = TRUE) {
  output_def <- c()

  args <- as.list(environment())
  for (a in names(args)) {
    if (a != "output_def") {
      if (args[[a]]) {
        output_def <- c(output_def, a)
      }
    }
  }
  return(output_def)
}

#' Let inits be randomly drawn using the heuristics by `stan`.
random_inits <- function(standat, model_type) {
  return(NULL)
}


#' Use inits informed by modeling assumptions, priors and/or data.
#'
#' Note that in contrast to using data to inform priors, using data to inform
#' inits does not result in double dipping, as it only defines a suitable
#' starting point for the stan sampler. Unless sampling is highly degenerate,
#' this does not bias the posterior. In contrast, starting closer to the target
#' region of the posterior can improve the sampling and avoid divergent
#' transitions etc.
#' 
#' The first part of ths function uses simple heuristic to get crude estimates
#' of important parameters. These are then used to define inits for these
#' parameters.
#'
#' @param standat A `list` with data that will be passed to stan (and which is
#'   partly also used by the inits here).
#' @param model_type Name of the type of model used.
#'
#' @return A function that can be called to obtain inits for one chain.
default_inits <- function(standat, model_type) {

  # get empirical estimates of cases and missingness
  known_by_reference_date <- rowSums(standat$reported_known)
  if (sum(standat$reported_unknown) > 0) {
    unknown_by_reference_date <- rep(0, standat$T)
    for (i in 1:standat$T) {
      # assume flat reporting delay
      unknown_fraction <- standat$reported_unknown[i] /
        (standat$D + 1)
      unknown_by_reference_date[max((i - standat$D), 1):i] <-
        unknown_by_reference_date[max((i - standat$D), 1):i] +
        unknown_fraction
    }
    all_by_reference_date <- unknown_by_reference_date + known_by_reference_date
    alpha_emp <- known_by_reference_date / all_by_reference_date
    alpha_emp <- mean(alpha_emp[1:7])
  } else {
    all_by_reference_date <- known_by_reference_date
    alpha_emp <- rep(1 - 1e-4, 7)
  }
  alpha_emp_logit <- qlogis(alpha_emp)
  # replace NAs with 50% probability (qlogis(0.5)=0)
  alpha_emp_logit <- ifelse(is.na(alpha_emp_logit), 0, alpha_emp_logit)

  # scale cases by assumed reporting proportion
  reporting_proportion <- sum(standat$latent_delay_dist)
  known_by_reference_date <- known_by_reference_date / reporting_proportion
  if (exists("unknown_by_reference_date")) {
    unknown_by_reference_date <-
      unknown_by_reference_date / reporting_proportion
  }
  all_by_reference_date <- all_by_reference_date / reporting_proportion

  # get crude estimate of baseline hazard, assume beta and eta are zero
  gamma_emp <- rowMeans(sapply(1:7, function(i) {
    standat$reported_known[i, ] /
      sum(standat$reported_known[i, ])
  }))
  gamma_emp_haz <- get_hazard_from_p(gamma_emp)
  gamma_emp_haz_reduced <- rep(0, standat$n_delays)
  for (i in 1:standat$n_delays) {
    gamma_emp_haz_reduced[i] <- mean(
      gamma_emp_haz[standat$delay_idx == i]
    )
  }

  # assume infections = symptom onsets shifted by mean incubation period
  medianInc <- weighted.median(
    0:(length(standat$latent_delay_dist) - 1),
    standat$latent_delay_dist
  )
  shift_backward <- function(x, shift) {
    c(x[(1 + shift):length(x)], rep(NA, times = shift))
  }
  I_length <- standat$L + standat$n_lambda_pre + standat$T
  start_t <- standat$L + standat$n_lambda_pre
  I_by_reference_date <- rep(NA, I_length)
  I_by_reference_date[
    (start_t + 1):(start_t + length(all_by_reference_date))
  ] <- all_by_reference_date
  firstI <- I_by_reference_date[start_t + 1]
  lastI <- I_by_reference_date[start_t + length(all_by_reference_date)]
  I_by_reference_date <- shift_backward(I_by_reference_date, medianInc)
  I_by_reference_date[1:(start_t - medianInc)] <- firstI
  I_by_reference_date[
    (start_t + length(all_by_reference_date) - medianInc):
    (start_t + length(all_by_reference_date))
  ] <- lastI
  last_unbiased <- length(I_by_reference_date) - standat$D
  I_by_reference_date <- c(
    I_by_reference_date[1:last_unbiased],
    rep(
      I_by_reference_date[last_unbiased],
      standat$D
    )
  )
  # back-calculate infections and R time series
  gt_mean <- weighted.mean(
    0:(length(standat$generation_time_dist) - 1),
    standat$generation_time_dist
  )
  gt_var <- weighted.mean(
    (0:(length(standat$generation_time_dist) - 1))^2,
    standat$generation_time_dist
  ) - gt_mean^2
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

    # use an empirical estimate of alpha and assume it remains
    # roughly unchanged over time
    inits[["alpha_logit_start"]] <- max(alpha_emp_logit, -10)
    inits[["alpha_logit_noise"]] <- rep(0, standat$T)
    inits[["alpha_logit_sd"]] <- max(
      standat$alpha_logit_sd_prior_mu +
        standat$alpha_logit_sd_prior_sd / 2, 1e-4
    )

    # reporting delay and effects
    inits[["gamma"]] <- gamma_emp_haz_reduced + 1e-4
    inits[["beta"]] <- rep(1e-4, standat$n_beta)
    inits[["beta_sd"]] <- 0.5
    inits[["eta"]] <- rep(1e-4, standat$n_eta)

    # overdispersion
    inits[["xi_negbinom"]] <- max(standat$xi_negbinom_prior_mu +
      standat$xi_negbinom_prior_sd / 2, 1e-4)

    # exponential smoothing
    inits[["ets_alpha"]] <- 1 - 1e-4
    inits[["ets_beta"]] <- 1 - 1e-4
    inits[["ets_phi"]] <- 1e-4
    
    # expected cases (if stepwise approach)
    if (model_type %in% c("impute_adjust", "adjust")) {
      inits[["lambda_log_start_values"]] <- c(
        log(max(standat$expected_cases_start, 1e-4)),
        rep(0, 1 + standat$ets_diff)
      )
      inits[["lambda_log_noise"]] <- rep(
        0, standat$T + standat$n_lambda_pre - 1 - standat$ets_diff
      )
      inits[["lambda_log_sd"]] <- max(
        standat$lambda_log_sd_prior_mu + standat$lambda_log_sd_prior_sd / 2,
        1e-4)
    }

    # R time series
    inits[["R_level_start"]] <- R_emp[standat$max_gen + 1]
    inits[["R_trend_start"]] <- 1e-4
    inits[["R_sd"]] <- 1
    inits[["R_noise"]] <- rep(
      0, standat$L + standat$n_lambda_pre + standat$T - standat$max_gen - 1
      )
    
    # infections
    # use initial expected cases as proxy for infections / expected infections
    inits[["iota_log_ar_start"]] <- log(max(I_by_reference_date[1], 1e-4))
    inits[["iota_log_ar_sd"]] <- 1
    inits[["iota_log_ar_noise"]] <- c(
      0, diff(log(pmax(I_by_reference_date, 1e-4)))[1:(standat$max_gen - 1)]
      )
    inits[["I"]] <- I_by_reference_date + 0.1
    return(inits)
  }

  return(init_function)
}

#' Ensure that stan model is up-to-date and potentially recompile
update_compiled_stanmodel <- function(model_def, code_dir = here::here("code"),
                                      force_recompile = FALSE) {
  n_models <- length(model_def$model_name)

  model_path_list <- lapply(1:n_models, function(i) {
    file.path(code_dir, model_def$model_folder, model_def$model_name[[i]])
  })
  include_paths_list <- lapply(1:n_models, function(i) {
    file.path(code_dir, model_def$model_folder) # identical
  })

  if (!model_def$profile) {
    stan_no_profiling_list <- lapply(1:n_models, function(i) {
      write_stan_files_no_profiling(
        model_path_list[[i]], include_paths_list[[i]]
      )
    })
    model_path_list <- lapply(1:n_models, function(i) {
      stan_no_profiling_list[[i]]$model
    })
    include_paths_list <- lapply(1:n_models, function(i) {
      stan_no_profiling_list[[i]]$include_paths
    })
  }

  cpp_options <- list()
  if (model_def$threads) {
    cpp_options[["stan_threads"]] <- TRUE
  }

  stanmodel_list <- lapply(1:n_models, function(i) {
    try(cmdstan_model(
      stan_file = model_path_list[[i]],
      include_paths = include_paths_list[[i]],
      dir = dirname(model_path_list[[i]]),
      cpp_options = cpp_options,
      force_recompile = force_recompile
    ))
  })

  model_def[["get_stan_model"]] <- lapply(1:n_models, function(i) {
    function() {
      return(stanmodel_list[[i]])
    }
  })

  return(model_def)
}
