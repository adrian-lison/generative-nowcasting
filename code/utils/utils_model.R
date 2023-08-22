#' Define details of the nowcasting model to be used.
#'
#' @param model_type The specific model to be used. One out of: ToDo
#' @param D The maximum assumed reporting delay
#' @param delay_resolution A `vector` describing the resolution at which the
#'   discrete reporting delay distribution should be modeled. The `vector`
#'   counts the number of days which are modeled together in one bin, starting
#'   with a delay of zero. The last entry of the `vector` must be Inf, to mark
#'   the maximum delay. For example, `c(1, rep(1, 14), 7, 7, Inf)` represents a
#'   distribution where delays of zero up to 14 days are modeled at a daily
#'   level, but the delays for 15-21 and for 22-28 are binned together as one
#'   week each. The sum of the `vector` (excluding Inf) must be equal to `D`.
#'   Default is `NULL`, i.e. daily resolution for all delays until the maximum
#'   delay.
#' @param delay_changepoint Changepoint model for reporting delay effects by
#'   date of reference. Either a random walk ("rw") or a piecewise linear
#'   ("segmented") model.
#' @param n_lambda_pre The number of days before `start_date` for which the
#'   reference date time series should already be modeled. The longer this time
#'   series is, the more report dates from `prep_data_missing` can be included
#'   in the model likelihood. At the same time, the dates modeled before
#'   `start_date` can be highly uncertain as they are only informed by cases
#'   with missing reference date.
#' @param overdispersion Should observations be modeled with (negative binomial)
#'   or without (Poisson) overdispersion?
#' @param ts_model What time series smoothing prior should be used to smooth
#'   latent time series such as the reproduction number over time? One of "none"
#'   (not recommendable unless you have lots of data), "ets" (exponential
#'   smoothing), or "sma" (simple MA model). Currently, "ets" is the default and
#'   only recommended option.
#' @param sma_window If using the "sma" smoothing prior, what should the window
#'   size for the MA component be?
#' @param ets_diff If using the "ets" prior, should the smoothing be applied to
#'   the actual time series or to the differenced time series? Default is FALSE.
#' @param ets_noncentered  if using the "ets" prior, should a non-centered
#'   parameterization be used for the noise/innovations component? Default is
#'   TRUE (recommended).
#' @param profile Should profiling information be collected during model
#'   fitting? Otherwise, a model with profiling statements removed is used, see
#'   `cmdstan_model_optional_profiling()`.
#' @param threads Should within-chain parallelization be supported using
#'   threading?
#' @param force_recompile Should recompilation be enforced?
#'
#' @return A `list` with the model type and details of the model definition.
define_model <- function(model_type = "base",
                         D = 28,
                         delay_resolution = NULL,
                         delay_changepoint = "rw",
                         n_lambda_pre = 1,
                         overdispersion = FALSE,
                         ts_model = "ets",
                         sma_window = 0,
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
    impute_forward = list("impute_forward.stan"),
    nowcast_imputed = list("nowcast_multiple_impute.stan"),
    nowcast_imputed_renewal = list("nowcast_multiple_impute_renewal.stan"),
    impute_and_nowcast = list("impute_nonparametric.stan",
                              "nowcast_multiple_impute.stan"),
    impute_forward_and_nowcast = list("impute_forward.stan",
                                      "nowcast_multiple_impute.stan"),
    impute_and_nowcast_renewal = list("impute_nonparametric.stan",
                                      "nowcast_multiple_impute_renewal.stan"),
    impute_forward_and_nowcast_renewal = list("impute_forward.stan",
                                              "nowcast_multiple_impute_renewal.stan"),
    impute_parametric_and_nowcast = list("impute_parametric_nb.stan",
                                         "nowcast_multiple_impute.stan"),
    renewal_direct = list("renewal.stan"),
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
    # last delay resolution must be Inf
    stopifnot(delay_resolution[length(delay_resolution)] == Inf)

    # remove last delay resolution (we do not model it)
    delay_resolution <- delay_resolution[-length(delay_resolution)]

    # ensure consistency with D
    stopifnot(sum(delay_resolution) == D)

    delay_idx <- rep(seq(1, length(delay_resolution)), times = delay_resolution)
    # add final elements for D and beyond
    delay_idx <- c(delay_idx, max(delay_idx) + 1, max(delay_idx) + 2)
  }
  model_def[["delay_idx"]] <- delay_idx

  model_def[["force_recompile"]] <- force_recompile
  model_def[["threads"]] <- threads
  model_def[["profile"]] <- profile

  model_def <- update_compiled_stanmodel(
    model_def, force_recompile = force_recompile
  )

  return(model_def)
}
