#' Define details of the nowcasting model to be used.
#' 
#' @details The following models are available:
#' - renewal_direct (renewal model for direct R estimation)                       
#' - impute_adjust (joint model for imputation and truncation adjustment)
#' - impute_adjust_latent (joint model for imputation and truncation adjustment)
#' - impute_adjust_renewal (joint model for imputation, truncation adjustment 
#' and R estimation)
#' - impute (imputation model for missing reference dates)
#' - impute_parametric (imputation model for missing reference dates, using a 
#' parametric instead of a non-parametric delay distribution model)
#' - impute_independent (imputation model for missing reference dates, using 
#' independent imputation, no backward delay model)
#' - adjust (model for truncation adjustment, supports multiple imputed data
#' as input)
#' - adjust_renewal (joint model for truncation adjustment and R estimation, 
#' supports multiple imputed data as input)
#' - impute_then_adjust (stepwise: imputation using model "impute", truncation 
#' adjustment using model "adjust")
#' - impute_then_adjust_renewal (stepwise: imputation using model "impute", 
#' truncation adjustment and R estimation using joint model "adjust_renewal")
#' - impute_independent_then_adjust (stepwise: imputation using model 
#' "impute_independent", truncation adjustment using model "adjust")
#' - impute_independent_then_adjust_renewal (stepwise: imputation using model 
#' "impute_independent", truncation adjustment and R estimation using joint 
#' model "adjust_renewal")
#' - impute_parametric_then_adjust (stepwise: imputation using model 
#' "impute_parametric", truncation adjustment using model "adjust")
#' 
#' @details Models of the type "adjust" without a renewal component will also
#' estimate R in an additional step, once using the model "renewal_direct" and 
#' once using EpiEstim. Uncertainty from the previous steps is accounted for
#' via resampling.
#'
#' @param model_type The specific model to be used. See details for the
#'   different model types available.
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
#' @param model_folder In which subfolder of the "code" folder is the model
#'   stored?
#'
#' @return A `list` with the model type and details of the model definition.
define_model <- function(model_type = "impute_adjust",
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
                         force_recompile = FALSE,
                         model_folder = "models") {
  model_def <- list()

  model_def[["model_type"]] <- model_type
  model_def[["model_name"]] <- switch(model_type,
    renewal_direct = list("renewal.stan"),                                  
    impute_adjust = list("impute_adjust.stan"),
    impute_adjust_latent = list("impute_adjust_latent.stan"),
    impute_adjust_renewal = list("impute_adjust_renewal.stan"),
    impute = list("impute_backward.stan"),
    impute_parametric = list("impute_backward_negbin.stan"),
    impute_independent = list("impute_independent.stan"),
    adjust = list("adjust.stan"),
    adjust_renewal = list("adjust_renewal.stan"),
    impute_then_adjust = list("impute_backward.stan",
                              "adjust.stan"),
    impute_then_adjust_renewal = list("impute_backward.stan",
                                      "adjust_renewal.stan"),
    impute_independent_then_adjust = list("impute_independent.stan",
                                          "adjust.stan"),
    impute_independent_then_adjust_renewal = list("impute_independent.stan", 
                                                  "adjust_renewal.stan"),
    impute_parametric_then_adjust = list("impute_backward_negbin.stan",
                                         "adjust.stan"),
    stop("Unknown model type.")
  )

  model_def[["model_folder"]] <- model_folder

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
