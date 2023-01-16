source(here::here("code", "utils", "utils_standata.R"))
source(here::here("code", "utils", "utils_model.R"))
source(here::here("code", "utils", "utils_priors.R"))
source(here::here("code", "utils", "utils_fit.R"))
source(here::here("code", "utils", "utils_job.R"))

#' Gather all relevant nowcast specifications, fit nowcasting model, and return results and diagnostics
make_nowcast <- function(data_def, model_def, prior_def, sampling_def, output_def, index = NULL) {
  ## Job management ----
  n_datasets <- length(data_def[["ll_data"]])
  n_dates <- length(data_def[["now"]])
  dataset_index <- 1
  date_index <- 1
  if ((n_datasets > 1) || (n_dates > 1)) {
    if (is.null(index)) {
      stop("Multiple time windows and/or datasets supplied. Please provide index")
    }
    dataset_index <- 1 + ((index - 1) %/% n_dates)
    date_index <- 1 + ((index - 1) %% n_dates)
  }
  
  ## Select data and priors for right dataset and date
  data_def$ll_data <- data_def$ll_data[[dataset_index]]
  data_def$now <- data_def$now[date_index]
  data_def$start_date <- data_def$start_date[date_index]
  
  if (is.list(prior_def$reporting_proportion) && length(prior_def$reporting_proportion)>1) {prior_def$reporting_proportion <- prior_def$reporting_proportion[[date_index]]}
  if (is.list(prior_def$latent_delay_dist) && length(prior_def$latent_delay_dist)>1) {prior_def$latent_delay_dist <- prior_def$latent_delay_dist[[date_index]]}
  if (is.list(prior_def$generation_time_dist) && length(prior_def$generation_time_dist)>1) {prior_def$generation_time_dist <- prior_def$generation_time_dist[[date_index]]}
  
  print(paste("Selected time window",
              data_def[["start_date"]], "to", data_def[["now"]],
              "from dataset", dataset_index))

  ## Stan data ----
  stan_prep <- get_stanprep(data_def, model_def, prior_def, sampling_def)

  ## Inits ----
  inits <- sampling_def[["inits_function"]](stan_prep$stan_data_list, model_def[["model_type"]])

  ## Model fitting ----
  nowcast_result <- list()

  nowcast_result[["model_type"]] <- model_def$model_type
  if (!is.null(index)) nowcast_result[["job_id"]] <- index
  nowcast_result[["dataset_index"]] <- dataset_index
  nowcast_result[["date_index"]] <- date_index
  nowcast_result[["stan_data_list"]] <- stan_prep$stan_data_list
  nowcast_result[["stan_prep_info"]] <- stan_prep$additional_info
  nowcast_result[["priors"]] <- prior_def
  # make sure model is compiled and up-to-date
  model_def <- update_compiled_stanmodel(model_def)
  stanmodel <- model_def[["get_stan_model"]][[1]]()
  print(paste("Running model executable:",stanmodel$exe_file()))
  fitted_model <- NULL
  
  try(fitted_model <- stanmodel$sample(
    data = stan_prep$stan_data_list,
    iter_warmup = sampling_def[["iter_warmup"]],
    iter_sampling = sampling_def[["iter_sampling"]],
    adapt_delta = sampling_def[["adapt_delta"]],
    step_size = sampling_def[["step_size"]],
    max_treedepth = sampling_def[["max_treedepth"]],
    chains = sampling_def[["chains"]],
    parallel_chains = sampling_def[["parallel_chains"]],
    threads_per_chain = sampling_def[["threads_per_chain"]],
    seed = sampling_def[["seed"]],
    refresh = sampling_def[["refresh"]],
    init = inits,
    show_messages = sampling_def[["show_messages"]]
  ))

  ## Direct postprocessing ----
  # envoke all lazy read methods to ensure data is loaded
  try(fitted_model$draws(), silent = TRUE)
  try(fitted_model$sampler_diagnostics(), silent = TRUE)
  try(fitted_model$init(), silent = TRUE)

  if (!is.null(fitted_model)) {
    nowcast_result[["fit"]] <- fitted_model
    if (!all(fitted_model$return_codes() == 0)) {
      warning("Fitting was not succesful. Returning model object for diagnostics.")
    }
  }
  # temporary save
  try({
    if (sampling_def[["temp_save"]]) {
      temp_rds_file <- tempfile(fileext = ".rds")
      fitted_model$save_object(file = temp_rds_file)
      print(paste("Model temporarily saved at", temp_rds_file))

      temp_rds_file <- here::here("data", "fitted_models", "last_fitted.rds")
      fitted_model$save_object(file = temp_rds_file)
      print(paste("Model temporarily saved at", temp_rds_file))
    }
  })

  ## Summary of results ----
  try(nowcast_result[["summary"]] <- summarize_fit(
    fitted_model, output_def, model_def[["model_type"]], data_def[["start_date"]],
    data_def[["now"]], stan_prep$stan_data_list))

  ## Diagnostics ----
  nowcast_result[["diagnostics"]] <- NULL
  nowcast_result[["diagnostic_summary"]] <- NULL
  try({
    nowcast_result[["diagnostics"]] <- fitted_model$cmdstan_diagnose()
    nowcast_result[["diagnostic_summary"]] <- fitted_model$diagnostic_summary()
  })

  nowcast_result[["profiles"]] <- NULL
  try(nowcast_result[["profiles"]] <- fitted_model$profiles(), silent = TRUE)

  nowcast_result[["output"]] <- NULL
  if (!is.null(fitted_model)) {
    try({
      nowcast_result[["output"]] <- fitted_model$output()
      if (length(fitted_model$output_files(include_failed = FALSE)) < fitted_model$num_chains()) {
        print("Not all chains have produced an output file.")
      }
    })
  }
  
  ## Output filtering ----
  # remove elements not in output_def
  for (e in names(nowcast_result)) {
    if (!(e %in% output_def)) {
      nowcast_result[e] <- NULL
    }
  }
  
  return(nowcast_result)
}

#' Impute missing data via backward delays, then fit a nowcasting model with multiple imputation approximation
fit_impute_and_nowcast <- function(data_def, model_def, prior_def, sampling_def, output_def, index = NULL) {
  ## Imputation ----
  data_def_impute <- define_data(ll_data = data_def$ll_data,
                                 now = data_def$now,
                                 start_date = data_def$start_date - model_def$D, # correcting for cutoff D
                                 holidays = data_def$holidays,
                                 imputed_posterior = NULL)
  
  model_def_impute <- define_model(
    model_type = "impute",
    D = model_def$D,
    delay_resolution = model_def$delay_resolution,
    delay_changepoint = model_def$delay_changepoint,
    n_lambda_pre = 0,
    overdispersion = model_def$delay_resolution,
    ts_model = model_def$ts_model,
    sma_window = model_def$sma_window,
    ets_diff = model_def$ets_diff,
    ets_noncentered = model_def$ets_noncentered,
    profile = model_def$profile,
    threads = model_def$threads,
    force_recompile = model_def$force_recompile)
  prior_def_impute <- define_priors(
    model_def_impute
  )
  
  output_def_impute <- define_output(
    fit = FALSE,
    stan_data_list = FALSE,
    stan_prep_info = FALSE,
    priors = FALSE,
    delays = FALSE
  )
  
  sampling_def_impute <- sampling_def
  sampling_def_impute$inits_function <- random_inits
  
  print("Fitting multiple imputation model.")
  imputed <- make_nowcast(data_def_impute, model_def_impute, prior_def_impute, sampling_def_impute, output_def, index)
  
  ## Nowcasting ----
  data_def_nowcast <- define_data(ll_data = data_def$ll_data,
                                  now = data_def$now,
                                  start_date = data_def$start_date,
                                  holidays = data_def$holidays,
                                  imputed_posterior = imputed$summary$imputed_posterior)
  
  nowcast_modeltype <- switch(model_def$model_type,
                              impute_and_nowcast = "nowcast_imputed",
                              impute_and_nowcast_renewal = "nowcast_imputed_renewal",
                              stop("Unknown model type."))
  
  model_def_nowcast <- define_model(
    model_type = nowcast_modeltype,
    D = model_def$D,
    delay_resolution = model_def$delay_resolution,
    delay_changepoint = model_def$delay_changepoint,
    ts_model = model_def$ts_model,
    sma_window = model_def$sma_window,
    n_lambda_pre = 0,
    overdispersion = model_def$overdispersion,
    ets_diff = model_def$ets_diff,
    ets_noncentered = model_def$ets_noncentered,
    profile = model_def$profile,
    threads = model_def$threads,
    force_recompile = model_def$force_recompile)
  
  print("Fitting nowcasting model.")
  nowcast <- make_nowcast(data_def_nowcast, model_def_nowcast, prior_def, sampling_def, output_def, index)
  
  nowcast$summary$mean_delay_backward <- imputed$summary$mean_delay_backward
  
  return(nowcast)
}
