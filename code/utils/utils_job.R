# Utils for nowcasting jobs ----
define_job <- function(job_name,
                       data_def,
                       model_def,
                       prior_def,
                       sampling_def = define_sampling(),
                       output_def = define_output()) {
  return(list(
    job_name = job_name,
    data_def = data_def,
    model_def = model_def,
    prior_def = prior_def,
    sampling_def = sampling_def,
    output_def = output_def
  ))
}

make_nowcast_from_job <- function(job, index = NULL) {
  if (nowcast_job$model_def$model_type %in% c("impute_and_nowcast", "impute_and_nowcast_renewal", "impute_parametric_and_nowcast")) {
    return(fit_impute_and_nowcast(job$data_def, job$model_def, job$prior_def, job$sampling_def, job$output_def, index = index))
  } else {
    return(make_nowcast(job$data_def, job$model_def, job$prior_def, job$sampling_def, job$output_def, index = index))
  }
}
