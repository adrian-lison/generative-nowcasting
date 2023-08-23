#' Define the data to be used for nowcasting
#'
#' @param ll_data Line list data. A `list` if multiple data sets (e.g. different
#'   simulation runs) should be fitted.
#' @param now Date at which the nowcast is made (only information to that date
#'   is available/observed). A `vector` if multiple dates should be fitted.
#' @param start_date Earliest symptom onset date to include in the nowcast.
#'   Together with `now`, this defines the data window. A `vector` if multiple
#'   dates should be fitted.
#' @param holidays A `vector` containing dates with (national) holidays which
#'   will be treated as Sundays by the reporting effect model.
#' @param imputed_posterior When using a stepwise approach with an imputation
#'   step, this is an array with imputed samples of missing dates, see
#'   `get_imputed_posterior()`.
#'
#' @return A `list` with the data.
define_data <- function(ll_data,
                        now,
                        start_date,
                        holidays = NULL,
                        imputed_posterior = NULL) {
  stopifnot(length(now) == length(start_date))
  data_def <- list()
  if (!(class(ll_data) == "list")) {
    ll_data <- list(ll_data)
  }
  data_def[["ll_data"]] <- ll_data
  data_def[["now"]] <- now
  data_def[["start_date"]] <- start_date
  data_def[["holidays"]] <- holidays
  if (!is.null(imputed_posterior)) {
    data_def[["imputed_posterior"]] <- imputed_posterior
  }
  return(data_def)
}

#' Get a preprocessed data list for stan, including priors and model settings
#'
#' @param data_def The data to be used for nowcasting, see `define_data()`.
#' @param model_def The model to be used for nowcasting, see `define_model()`.
#' @param prior_def The priors to be used for nowcasting, see `define_priors()`.
#' @param sampling_def The sampling parameters to be used by stan, 
#' see `define_sampling()`.
#'
#' @return A `list` with the standat (`stan_data_list`) and the meta data 
#' (`additional_info`).
get_stanprep <- function(data_def, model_def, prior_def, sampling_def) {
  set.seed(sampling_def[["seed"]])

  # Make data list for stan
  prep_data_complete <- data_def[["ll_data"]] %>%
    filter(event2_date <= data_def[["now"]], !is.na(event1_date))
  prep_data_missing <- data_def[["ll_data"]] %>%
    filter(event2_date <= data_def[["now"]], is.na(event1_date))
  stan_prep <- prepare_data_list(
    prep_data_complete = prep_data_complete,
    prep_data_missing = prep_data_missing,
    now = data_def[["now"]],
    start_date = data_def[["start_date"]],
    holidays = data_def[["holidays"]],
    D = model_def[["D"]],
    delay_idx = model_def[["delay_idx"]],
    delay_changepoint = model_def[["delay_changepoint"]],
    n_lambda_pre = model_def[["n_lambda_pre"]],
    latent_delay_dist = prior_def[["latent_delay_dist"]],
    generation_time_dist = prior_def[["generation_time_dist"]],
    reporting_proportion = prior_def[["reporting_proportion"]]
  )

  ## add further structural model details
  stan_prep$stan_data_list[["overdispersion"]] <- model_def[["overdispersion"]]
  stan_prep$stan_data_list[["ts_model"]] <- model_def[["ts_model"]]
  stan_prep$stan_data_list[["sma_window"]] <- model_def[["sma_window"]]
  stan_prep$stan_data_list[["ets_diff"]] <- model_def[["ets_diff"]]
  stan_prep$stan_data_list[["ets_noncentered"]] <- model_def[["ets_noncentered"]]

  # add information about multiple imputations
  if (model_def$model_type %in% c("adjust", "adjust_renewal")) {
    if (!("imputed_posterior" %in% names(data_def))) {
      stop("Multiple imputed data missing.")
    } else {
      stan_prep$stan_data_list[["n_imputations"]] <-
        dim(data_def$imputed_posterior)[1]

      stan_prep$stan_data_list[["reported_unknown_imputed"]] <-
        data_def$imputed_posterior
    }
  }

  # add priors
  stan_prep$stan_data_list[["ets_fixed"]] <- prior_def[["ets_fixed"]]
  priors <- prior_def[["get_priors"]](stan_prep$stan_data_list)
  stan_prep$additional_info[["priors"]] <- priors
  stan_prep$stan_data_list <- append(
    stan_prep$stan_data_list,
    purrr::flatten(priors)
  )

  return(stan_prep)
}


#' Prepprocess data for nowcasting
#'
#' @param prep_data_complete Line list data of cases with known reference dates
#' @param prep_data_missing Line list data of cases with missing reference dates
#' @param now Date at which the nowcast is made (only information to that date
#'   is available/observed)
#' @param start_date Earliest symptom onset date to include in the nowcast.
#'   Together with `now`, this defines the data window.
#' @param holidays A `vector` containing dates with (national) holidays which
#'   will be treated as Sundays by the reporting effect model.
#' @param D Maximum assumed reporting delay
#' @param delay_idx Indexes for sparse delay specification
#' @param delay_changepoint Changepoint model for reporting delay effects by
#'   date of reference. Either a random walk ("rw") or a piecewise linear
#'   ("segmented") model.
#' @param changepoint_interval At what interval (in days) should the reporting
#'   delay effects be allowed to change?
#' @param n_lambda_pre The number of days before `start_date` for which the
#'   reference date time series should already be modeled. The longer this time
#'   series is, the more report dates from `prep_data_missing` can be included
#'   in the model likelihood. At the same time, the dates modeled before
#'   `start_date` can be highly uncertain as they are only informed by cases
#'   with missing reference date.
#' @param reporting_proportion Proportion of infections which become reported
#'   cases. 1 means full reporting (100% of infections are reported).
#' @param latent_delay_dist A `vector` representing a discrete distribution of
#'   the delay between latent (e.g. infections) and observed (e.g. symptom
#'   onsets) events, e.g. the incubation period.
#' @param generation_time_dist A `vector` representing a discrete generation
#'   time distribution.
#'
#' @return A `list` with the preprocessed data
prepare_data_list <- function(prep_data_complete,
                              prep_data_missing,
                              now,
                              start_date,
                              holidays = NULL,
                              D,
                              delay_idx,
                              delay_changepoint = "segmented",
                              changepoint_interval = 7,
                              n_lambda_pre,
                              latent_delay_dist,
                              generation_time_dist,
                              reporting_proportion = 1) {
  stan_data_list <- list() # directly used in the Stan model
  additional_info <- list() # used by other components / info for the user

  ### Durations
  all_dates <- seq(start_date, now, by = "1 day")
  # time series including the days of the modeling phase
  tswp <- seq(start_date - n_lambda_pre, now, by = "1 day")
  T <- length(all_dates) - 1
  pre_diff <- D - n_lambda_pre

  stan_data_list[["T"]] <- T + 1
  stan_data_list[["D"]] <- D
  stan_data_list[["n_lambda_pre"]] <- n_lambda_pre
  additional_info[["all_dates"]] <- all_dates
  stan_data_list[["n_delays"]] <- max(delay_idx) - 2
  stan_data_list[["delay_idx"]] <- delay_idx

  ### Reporting triangle
  reporting_data <- prepare_reporting_data(
    prep_data_complete, prep_data_missing, now, start_date, D
  )
  stan_data_list[["cases"]] <- reporting_data[["cases"]]
  stan_data_list[["reported_known"]] <- reporting_data[["reported_known"]]
  stan_data_list[["reported_unknown"]] <- reporting_data[["reported_unknown"]]
  if (pre_diff > 0) {
    stan_data_list[["emp_alpha"]] <- reporting_data[["emp_alpha"]][1:pre_diff]
  }

  stan_data_list[["emp_forward_delay_dist"]] <-
    reporting_data[["emp_forward_delay_dist"]]

  stan_data_list[["emp_backward_delay_dist"]] <-
    reporting_data[["emp_backward_delay_dist"]]

  # estimate initial expected cases based on
  # reported known + unknown (assuming a flat delay distribution)
  expected_cases_start <- sum(reporting_data[["reported_known"]][1, ]) +
    mean(reporting_data[["reported_unknown"]][1:D])
  stan_data_list[["expected_cases_start"]] <- expected_cases_start

  ### Date of occurrence effects on reporting hazard
  # Design matrix for delay effects by occurrence date
  if (delay_changepoint == "rw") {
    Z <- get_design_changepoint(
      start_date = start_date,
      now = now,
      n_lambda_pre = n_lambda_pre,
      change_interval = changepoint_interval,
      segmented = FALSE,
      scale = FALSE
    )
  } else if (delay_changepoint == "segmented") {
    Z <- get_design_changepoint(
      start_date = start_date,
      now = now,
      n_lambda_pre = n_lambda_pre,
      change_interval = changepoint_interval,
      segmented = TRUE,
      scale = TRUE
    )
  } else {
    stop("No valid delay changepoint model provided.", call. = F)
  }
  stan_data_list[["Z"]] <- Z
  stan_data_list[["n_beta"]] <- dim(Z)[2]

  ### Date of report effects on reporting hazard
  # Make designmatrix for weekdays
  # Make a factor for each day of the week
  wdays <- levels(lubridate::wday(tswp, label = TRUE))[-1]
  # Create the extension of the W matrix
  W <- array(NA,
    dim = c(length(tswp), D, length(wdays)),
    dimnames = list(
      as.character(tswp),
      paste("delay", 0:(D - 1), sep = ""), wdays
    )
  )
  # Replace holidays with sunday indicator
  get_wdays <- function(nc_dates_t, D) {
    wdays <- lubridate::wday(nc_dates_t + 0:(D - 1), label = TRUE)
    if (!is.null(holidays)) {
      wdays[(nc_dates_t + 0:(D - 1)) %in% holidays] <- "Sun"
    }
    return(wdays)
  }
  # Assign indicator to all times and lags
  for (t in seq_len(length(tswp))) {
    for (w in seq_len(length(wdays))) {
      #
      W[t, , w] <- as.numeric(get_wdays(tswp[t], D) == wdays[w])
    }
  }
  # Reduce design matrices
  if (!all(delay_idx == seq(1, D + 2))) {
    W <- reduce_design_matrix(W, delay_idx)
  }
  stan_data_list[["W"]] <- W
  stan_data_list[["n_eta"]] <- dim(W)[3]

  ### Delay distributions
  # Incubation period distribution (scaled by reporting proportion)
  latent_delay_dist <- latent_delay_dist * reporting_proportion
  stan_data_list[["L"]] <- length(latent_delay_dist) - 1 # includes zero
  stan_data_list[["latent_delay_dist"]] <- latent_delay_dist

  # Generation time distribution
  stan_data_list[["max_gen"]] <- length(generation_time_dist)
  stan_data_list[["generation_time_dist"]] <- generation_time_dist

  return(list(
    stan_data_list = stan_data_list,
    additional_info = additional_info
  ))
}

#' Get reporting triangles and empirical estimates of delays and missingness
#' 
#' @param prep_data_complete Line list data of cases with known reference dates
#' @param prep_data_missing Line list data of cases with missing reference dates
#' @param now Date at which the nowcast is made (only information to that date
#'   is available/observed)
#' @param start_date Earliest symptom onset date to include in the nowcast.
#'   Together with `now`, this defines the data window.
#' @param D Maximum assumed reporting delay
prepare_reporting_data <- function(prep_data_complete, prep_data_missing,
                                   now, start_date, D) {
  if (!("n" %in% names(prep_data_complete))) {
    prep_data_complete <- prep_data_complete %>% mutate(n = 1)
  }
  if (!("n" %in% names(prep_data_missing))) {
    prep_data_missing <- prep_data_missing %>% mutate(n = 1)
  }

  # Cases by date of report
  cases <- bind_rows(prep_data_complete, prep_data_missing) %>%
    group_by(event2_date) %>%
    dplyr::summarize(n = sum(n), .groups = "drop") %>%
    filter(event2_date >= start_date, event2_date <= now) %>%
    complete(
      event2_date = seq.Date(start_date, now, by = "1 day"),
      fill = list(n = 0)
    ) %>%
    pull(n)

  # Cases with observed type 1 event
  reporting_triangle <- prep_data_complete %>%
    filter(event1_date >= start_date, event2_date <= now) %>%
    mutate(delay = as.numeric(event2_date - event1_date))

  # Currently, we exclude observations with negative delay. If alternative
  # behavior is wanted (e.g. trim to zero, or set missing), this should already
  # be done before calling this function. Observations with delay longer than
  # the specified maximum delay are truncated.
  n_negative_delay <- reporting_triangle %>%
    ungroup() %>%
    filter(delay < 0) %>%
    summarize(n = sum(n), .groups = "drop") %>%
    pull(n)
  n_maximum_delay <- reporting_triangle %>%
    ungroup() %>%
    filter(delay > D) %>%
    summarize(n = sum(n), .groups = "drop") %>%
    pull(n)
  if (n_negative_delay > 0) {
    warning(paste("Removed ",
      n_negative_delay,
      " case(s) with negative delay.",
      sep = ""
    ))
  }
  if (n_maximum_delay > 0) {
    warning(paste("Removed ",
      n_maximum_delay,
      " cases with delay > ",
      D, ".",
      sep = ""
    ))
  }

  reporting_triangle %<>%
    filter(delay >= 0) %>% # remove negative delays
    filter(delay <= D) # remove delays longer than D

  # Compute reporting triangle
  reporting_triangle %<>%
    group_by(event1_date, delay) %>%
    dplyr::summarize(n = sum(n), .groups = "drop") %>%
    complete(delay = 0:D, event1_date, fill = list(n = 0)) %>%
    arrange(delay, event1_date) %>%
    pivot_wider(names_from = delay, values_from = n) %>%
    ungroup() %>%
    complete(event1_date = seq.Date(start_date, now, by = "1 day")) %>%
    mutate(across(-event1_date, replace_na, replace = 0))

  # Convert to matrix
  reporting_matrix <- as.matrix(reporting_triangle %>% select(-event1_date))
  dimnames(reporting_matrix)[[1]] <- as.character(reporting_triangle$event1_date)

  # Cases with unobserved type 1 event
  unobserved <- prep_data_missing %>%
    group_by(event2_date) %>%
    dplyr::summarize(n = sum(n), .groups = "drop") %>%
    filter(event2_date >= start_date, event2_date <= now) %>%
    complete(
      event2_date = seq.Date(start_date, now, by = "1 day"),
      fill = list(n = 0)
    )

  unobserved_matrix <- unobserved %>% pull(n)

  # Empirical estimates for seeding etc
  # rough estimate of forward delay, does not account for right-truncation
  emp_forward_delay_dist <- prep_data_complete %>%
    filter(event1_date >= start_date, event1_date <= now) %>%
    get_empirical_delay(D)

  emp_backward_delay_dist <- prep_data_complete %>%
    filter(event2_date >= start_date, event2_date <= now) %>%
    get_empirical_backward_delay(D)

  emp_alpha <- unobserved %>%
    inner_join(emp_backward_delay_dist, by = "event2_date") %>%
    mutate(a = n * p, event1_date = event2_date - delay) %>%
    group_by(event1_date) %>%
    summarize(unobserved = sum(a)) %>%
    inner_join(
      reporting_triangle %>%
        transmute(event1_date, observed = rowSums(across(-event1_date))),
      by = "event1_date"
    ) %>%
    mutate(alpha = observed / (observed + unobserved)) %>%
    ungroup() %>%
    mutate(alpha = ifelse(is.nan(alpha), NA, alpha)) %>%
    fill(alpha, .direction = "downup")

  return(list(
    cases = cases,
    reported_known = reporting_matrix,
    reported_unknown = unobserved_matrix,
    emp_forward_delay_dist = emp_forward_delay_dist$p,
    emp_backward_delay_dist = emp_backward_delay_dist$p,
    emp_alpha = emp_alpha$alpha
  ))
}

## Helper functions ----
#' Get overall empirical delay from a line list, no smoothing and
#' no adjustment for right-truncation
get_empirical_delay <- function(df, D) {
  df %>%
    mutate(delay = as.numeric(event2_date - event1_date)) %>%
    filter(delay >= 0) %>% # remove negative delays
    mutate(delay = ifelse(delay > D, D, delay)) %>% # truncate delays longer than D
    count(delay) %>%
    complete(delay = 0:D, fill = list(n = 0)) %>%
    transmute(delay, p = n / sum(n)) %>%
    return()
}

#' Get daily empirical backward delay from a line list, no smoothing
get_empirical_backward_delay <- function(df, D) {
  df %>%
    mutate(delay = as.numeric(event2_date - event1_date)) %>%
    filter(delay >= 0) %>% # remove negative delays
    mutate(delay = ifelse(delay > D, D, delay)) %>% # truncate delays longer than D
    count(event2_date, delay) %>%
    complete(event2_date, delay = 0:D, fill = list(n = 0)) %>%
    group_by(event2_date) %>%
    transmute(delay, p = n / sum(n)) %>%
    return()
}

#' Compute backward delay from a time series and forward delay distribution
compute_backward_delay <- function(time_series, forward_delay_dist) {
  n_delays <- dim(forward_delay_dist)[2]
  bd <- matrix(0, nrow = length(time_series) + n_delays - 1, ncol = n_delays)
  for (i in 1:length(time_series)) {
    for (d in 1:n_delays) {
      bd[i + d - 1, d] <- bd[i + d - 1, d] + time_series[i] * forward_delay_dist[i, d]
    }
  }
  for (i in 1:(length(time_series) + n_delays - 1)) {
    if (sum(bd[i, ]) > 0) bd[i, ] <- bd[i, ] / sum(bd[i, ])
  }
  return(bd)
}

#' Get design matrix for changepoint model
get_design_changepoint <- function(start_date, now, n_lambda_pre,
                                   segmented = FALSE, change_interval = 7,
                                   scale = TRUE) {
  all_dates <- seq(start_date, now, by = "1 day")
  # time series including the days of the modeling phase
  tswp <- seq(start_date - n_lambda_pre, now, by = "1 day")
  T <- length(all_dates) - 1
  ddChangepoint <- rev(seq.Date(now, start_date + 1,
    by = str_glue("-{change_interval} days")
  ))
  Z_idx <- rev(rep(length(ddChangepoint):1, each = change_interval)[1:(1 + T)])
  Z <- array(NA,
    dim = c(T + 1 + n_lambda_pre, length(ddChangepoint)),
    dimnames = list(as.character(tswp), 1:length(ddChangepoint))
  )

  for (t_pre in 1:n_lambda_pre) {
    # we don't have a changepoint model before the inference phase, but fill
    # the space up with 0s to have the same dimensionality as the other matrices
    Z[t_pre, ] <- rep(0, length(ddChangepoint))
  }

  if (segmented) {
    # use segmented linear model
    for (t in 0:T) {
      covariates_ind_days <- c(0, as.numeric(all_dates <= (all_dates[t + 1]))[-1])
      # add up counts in each changepoint interval
      Z[n_lambda_pre + 1 + t, ] <- tapply(covariates_ind_days, Z_idx, sum)
    }
  } else {
    # use random walk
    for (t in 0:T) {
      covariates_ind_days <- c(TRUE, (all_dates <= (all_dates[t + 1]))[-1])
      # add up counts in each changepoint interval
      Z[n_lambda_pre + 1 + t, ] <- tapply(covariates_ind_days, Z_idx, any)
    }
  }

  if (all(Z[, 1] == 1)) {
    Z <- Z[, -1]
  }

  if (scale) {
    Z <- Z / (max(Z) * 0.68) # scale columns by maximum for better efficiency
  }
  return(Z)
}

#'  Reduction of design matrices for lower resolution intervals
reduce_design_matrix <- function(design_matrix, delay_idx) {
  n_time_points <- dim(design_matrix)[1]
  n_reduced_times <- max(delay_idx) - 2
  n_covariates <- dim(design_matrix)[3]
  reduced_matrix <- array(dim = c(n_time_points, n_reduced_times, n_covariates))

  for (d in 1:n_reduced_times) {
    delay_indices <- delay_idx[1:(length(delay_idx) - 2)] == d
    if (sum(delay_indices) == 1) {
      reduced_matrix[, d, ] <- design_matrix[, delay_indices, ]
    } else {
      reduced_matrix[, d, ] <- t(sapply(1:n_time_points, function(t) {
        colMeans(design_matrix[t, delay_indices, ])
      }))
    }
  }
  return(reduced_matrix)
}
