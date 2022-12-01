# Utils for selecting / defining the nowcasting model ----
define_model <- function(model_type = "base",
                         D = 28,
                         delay_resolution = NULL,
                         delay_changepoint = "rw",
                         n_lambda_pre = 1,
                         overdispersion = FALSE,
                         ts_model = "ets",
                         sma_window = 7,
                         ets_diff = FALSE,
                         ets_noncentered = TRUE,
                         profile = FALSE,
                         threads = FALSE,
                         force_recompile = FALSE) {
  model_def <- list()

  model_def[["model_type"]] <- model_type
  model_def[["model_name"]] <- switch(model_type,
    base = list("nowcast_impute.stan"),
    base_old = list("nowcast_impute.stan"),
    latent = list("nowcast_latent_impute.stan"),
    renewal = list("nowcast_renewal_impute.stan"),
    renewal_noncentered = list("nowcast_renewal_impute_noncentered.stan"),
    renewal_deterministic = list("nowcast_renewal_impute_deterministic.stan"),
    impute = list("impute_nonparametric.stan"),
    impute_parametric = list("impute_parametric_nb.stan"),
    nowcast_imputed = list("nowcast_multiple_impute.stan"),
    nowcast_imputed_renewal = list("nowcast_multiple_impute_renewal.stan"),
    impute_and_nowcast = list("impute_nonparametric.stan", "nowcast_multiple_impute.stan"),
    impute_and_nowcast_renewal = list("impute_nonparametric.stan", "nowcast_multiple_impute_renewal.stan"),
    impute_parametric_and_nowcast = list("impute_parametric_nb.stan", "nowcast_multiple_impute.stan"),
    stop("Unknown model type.")
  )
  if (model_type %in% c("base_old")) {
    model_def[["model_folder"]] <- "models_old"
  } else {
    model_def[["model_folder"]] <- "models"
  }

  model_def[["D"]] <- D
  model_def[["delay_resolution"]] <- delay_resolution
  model_def[["delay_changepoint"]] <- delay_changepoint
  model_def[["n_lambda_pre"]] <- n_lambda_pre
  model_def[["overdispersion"]] <- overdispersion
  if (ts_model %in% c("none", 0)) {
    model_def[["ts_model"]] <- 0 # no smoothing
  } else if (ts_model %in% c("ets", 1)) {
    model_def[["ts_model"]] <- 1 # exponential smoothing
  } else if (ts_model %in% c("sma", 2)) {
    model_def[["ts_model"]] <- 2 # simple MA model
  } else {
    stop("No valid time series model provided.", call. = F)
  }
  model_def[["sma_window"]] <- sma_window
  model_def[["ets_diff"]] <- ets_diff
  model_def[["ets_noncentered"]] <- ets_noncentered

  # delay index
  if (is.null(delay_resolution)) {
    delay_idx <- seq(1, D + 2)
  } else {
    stopifnot(delay_resolution[length(delay_resolution)] == Inf) # last delay resolution must be Inf
    delay_resolution <- delay_resolution[-length(delay_resolution)] # remove last delay resolution (we do not model it)
    stopifnot(sum(delay_resolution) == D) # ensure consistency with D
    delay_idx <- rep(seq(1, length(delay_resolution)), times = delay_resolution)
    delay_idx <- c(delay_idx, max(delay_idx) + 1, max(delay_idx) + 2) # add final elements for D and beyond
  }
  model_def[["delay_idx"]] <- delay_idx

  model_def[["force_recompile"]] <- force_recompile
  model_def[["threads"]] <- threads
  model_def[["profile"]] <- profile

  model_def <- update_compiled_stanmodel(model_def, force_recompile)

  return(model_def)
}
