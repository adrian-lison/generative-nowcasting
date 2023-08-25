#' Define a nowcasting job(array)
#' 
#' A nowcasting job contains all the relevant information (data, model, priors,
#' sampling parameters etc.) to fit a nowcasting model using `cmdstanr`. If 
#' using the same model over a number of dates or datasets with otherwise 
#' identical settings, a job array can be defined (this is achieved by 
#' providing lists of nowcast dates or linelist data to `data_def()`).
#'
#' @param job_name The name of the job(array).
#' @param data_def The data to be used for nowcasting, see `define_data()`.
#' @param model_def The model to be used for nowcasting, see `define_model()`.
#' @param prior_def The priors to be used for nowcasting, see `define_priors()`.
#' @param sampling_def The sampling parameters to be used by stan, 
#' see `define_sampling()`.
#' @param output_def The results to be returned from the fitted model, 
#' see `define_output()`.
#' @param index_by_date Should the job result files also be indexed by the 
#' nowcast date?
#' @param verbose Should a message about the job(array) be printed?
define_job <- function(job_name,
                       data_def,
                       model_def,
                       prior_def,
                       sampling_def = define_sampling(),
                       output_def = define_output(),
                       index_by_date = F,
                       verbose = T) {
  n_datasets <- length(data_def[["ll_data"]])
  n_dates <- length(data_def[["now"]])

  if (verbose) {
    print(paste0(
      "Creating nowcast job(array) ''",
      job_name, "'' of ",
      length(data_def$ll_data),
      " different datasets with ",
      length(data_def$now),
      " different dates each."
    ))
  }

  return(list(
    job_name = job_name,
    data_def = data_def,
    model_def = model_def,
    prior_def = prior_def,
    sampling_def = sampling_def,
    output_def = output_def,
    index_by_date = index_by_date,
    jobindex_mapping = 1:(n_datasets * n_dates)
  ))
}

#' Fit nowcasting model for nowcast job locally
#'
#' @param job A nowcast job, see `define_job()`.
#' @param index If part of a job array, the index of the job to be fitted.
#'
#' @return A `list` with nowcasting results, see `make_nowcast()` and 
#' `define_output()`.
make_nowcast_from_job <- function(job, index = NULL) {
  imputation_models <- c(
    "impute_then_adjust",
    "impute_independent_then_adjust",
    "impute_then_adjust_renewal",
    "impute_independent_then_adjust_renewal",
    "impute_parametric_then_adjust"
  )

  if (job$model_def$model_type %in% imputation_models) {
    return(fit_impute_and_nowcast(
      job$data_def,
      job$model_def,
      job$prior_def,
      job$sampling_def,
      job$output_def,
      index = index
    ))
  } else {
    return(make_nowcast(
      job$data_def,
      job$model_def,
      job$prior_def,
      job$sampling_def,
      job$output_def,
      index = index
    ))
  }
}
