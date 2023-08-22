#' Define simulation parameters
#'
#' @param initial_infections A vector as long as the maximum generation time
#'   with initial/seeding infections
#'
#' @param R A vector of the instantaneous effective reproduction number by date
#'   of infection
#'
#' @param generation_dist A vector of the generation time distribution. Must sum
#'   to one.
#'
#' @param incubation_dist A vector of the incubation time distribution. Must sum
#'   to one.
#'
#' @param hosp_prob A vector of the probability of hospitalization. by date of
#'   infection.
#'
#' @param symp_to_hosp_dist A vector of the time distribution from symptom onset
#'   to hospitalization. Must sum to one.
#'
#' @param onset_known_prob A vector of the probability of having the symptom
#'   onset recorded, by date of symptom onset.
#'
#' @param symp_to_rep_dist A vector of the time distribution from symptom onset
#'   to reporting. Must sum to one.
#'
#' @param hosp_to_rep_dist A vector of the time distribution from
#'   hospitalization to reporting. Must sum to one.
define_simulation <- function(
    initial_infections, R, generation_dist, incubation_dist, hosp_prob,
    symp_to_hosp_dist, onset_known_prob, symp_to_rep_dist, hosp_to_rep_dist) {
  simulation_def <- list()
  simulation_def[["initial_infections"]] <- initial_infections
  simulation_def[["R"]] <- R
  simulation_def[["generation_dist"]] <- generation_dist
  simulation_def[["incubation_dist"]] <- incubation_dist
  simulation_def[["hosp_prob"]] <- hosp_prob
  simulation_def[["symp_to_hosp_dist"]] <- symp_to_hosp_dist
  simulation_def[["onset_known_prob"]] <- onset_known_prob
  simulation_def[["symp_to_rep_dist"]] <- symp_to_rep_dist
  simulation_def[["hosp_to_rep_dist"]] <- hosp_to_rep_dist
  return(simulation_def)
}

#' Simulate ground truth linelist data based on a renewal process
simulate_linelist <- function(simulation_def) {
  infections <- simulate_infections(
    simulation_def$initial_infections,
    simulation_def$R,
    simulation_def$generation_dist
  )
  print(paste("Overall number of infections:", sum(infections[["I"]])))
  cases <- simulate_outcomes_cases(
    infections[["I"]],
    simulation_def$incubation_dist,
    simulation_def$hosp_prob,
    simulation_def$symp_to_hosp_dist
  )
  cases <- simulate_ascertainment_cases(
    cases,
    simulation_def$onset_known_prob,
    simulation_def$symp_to_rep_dist,
    simulation_def$hosp_to_rep_dist
  )
  process_summary <- summarize_linelist(cases, 1:length(infections[["I"]]))
  process_summary[["infections"]] <- infections[["I"]]
  process_summary[["infections_expected"]] <- infections[["I_expected"]]
  process_summary[["t_onset"]] <- 1:(length(infections[["I"]]) +
    ncol(simulation_def$incubation_dist) - 1)
  return(list(
    parameters = simulation_def,
    linelist = cases,
    process_summary = process_summary
  ))
}

#' Simulate infections using a renewal process
#' 
#' If an upper and lower number of infections per day are specified,
#' realizations of the scenario are sampled via rejection sampling.
#'
#' @param initial_infections A vector of initial/seeding infections that is as
#' long as the maximum generation time
#'
#' @param R A vector of the effective reproduction number over time. The length
#' of the vector determines the number of days infections are simulated.
#'
#' @param generation_dist A `matrix` (`(length(R)+maxGen)`x`maxGen`) of the
#' generation time distribution over time. In each row, the probabilities for
#' the discrete generation times from one to the maximum generation time are
#' given. Each row stands for a different day of primary infection (starting
#' with the initial infections, and then one for each day of `R`).
#' 
#' @param upper Maximum number of infections per day.
#' 
#' @param lower Minimum number of infections per day.
#'
#' @return A `list` with vectors for the expected infections and the realized
#' infections, from `1:length(R)` respectively.
#'
simulate_infections <- function(initial_infections, R, generation_dist,
                                upper = 25000, lower = 250) {
  maxGen <- length(generation_dist[1, ])
  assertthat::are_equal(length(initial_infections), maxGen)
  assertthat::are_equal(maxGen + length(R), dim(generation_dist)[1])

  max_inf <- as.numeric("Inf")
  while ((max_inf > upper) || (max_inf < lower)) {
    I_expected <- rep(0, length(R))
    I <- c(initial_infections, rep(0, length(R)))
    for (t in 1:length(R)) {
      # renewal equation
      I_expected[t] <- R[t] * sum(
        I[(maxGen + t - 1):(maxGen + t - maxGen)] *
          generation_dist[
            cbind((maxGen + t - 1):(maxGen + t - maxGen), 1:maxGen)
          ]
      )
      I[maxGen + t] <- rpois(1, I_expected[t])
    }
    max_inf <- max(I)
  }
  return(list(I_expected = I_expected, I = I[(maxGen + 1):length(I)]))
}

simulate_outcomes_cases <- function(
    infections, incubation_dist, hosp_prob, symp_to_hosp_dist) {
  infections_hosp <- sapply(
    1:length(infections),
    # sample how many infections are hospitalized, assuming probability of
    # hospitalization is time-dependent with respect to the time of infection
    function(t) rbinom(1, infections[t], hosp_prob[t])
  )
  onsets <- matrix(sapply(
    1:length(infections_hosp),
    # simulate incubation times
    function(t) rmultinom(1, infections_hosp[t], incubation_dist[t, ])
  ), ncol = 1)
  onsets <- cbind(
    # add infection times
    rep(1:length(infections_hosp),
      each = length(incubation_dist[1, ])
    ), onsets
  )
  onsets <- cbind(
    onsets[, 1],
    # compute onset times
    onsets[, 1] + rep(0:(length(incubation_dist[1, ]) - 1),
      times = length(infections_hosp)
    ), onsets[, 2]
  )
  colnames(onsets) <- c("infection_time", "onset_time", "onset_count")
  onsets <- as.data.frame(onsets)

  stopifnot(length(hosp_prob) >= max(onsets$infection_time))
  stopifnot(dim(symp_to_hosp_dist)[1] >= max(onsets$onset_time))

  individual_cases <- onsets %>%
    uncount(onset_count) %>%
    # we assume that the distribution of the time from symptom onset to
    # hospitalization is time-dependent w.r.t. the time of symptom onset
    mutate(hosp_time = sapply(
      onset_time,
      function(t) t + extraDistr::rcat(1, symp_to_hosp_dist[t, ])
    ))

  return(individual_cases)
}

#' Simulate the ascertainment / reporting process for each individual case.
#'
#' @param individual_cases Line list of individual cases with disease onset.
#' @param onset_known_prob The (time-varying) probability that the disease onset
#'   will be known (i.e. not missing) in the reporting data.
#' @param symp_to_rep_dist The (time-varying) delay distribution for the time
#'   from symptom onset to reporting.
#' @param hosp_to_rep_dist The (time-varying) delay distribution for the time
#'   from hospitalization to reporting.
#'
#' @return Line list of individual cases with disease onset and reporting data.
simulate_ascertainment_cases <- function(
    individual_cases, onset_known_prob,
    symp_to_rep_dist = NULL, hosp_to_rep_dist = NULL) {
  if (is.null(symp_to_rep_dist) & is.null(hosp_to_rep_dist)) {
    stop("Either delay from symptom onset or from hospitalization must be provided.")
  }
  stopifnot(length(onset_known_prob) >= max(individual_cases$onset_time))
  if (is.null(symp_to_rep_dist)) {
    delay_dist <- hosp_to_rep_dist
    stopifnot(dim(delay_dist)[1] >= max(individual_cases$hosp_time, na.rm = T))
  } else {
    delay_dist <- symp_to_rep_dist
    stopifnot(dim(delay_dist)[1] >= max(individual_cases$onset_time, na.rm = T))
  }
  individual_cases <- individual_cases %>%
    # we assume that the probability of the symptom onset being recorded is
    # time-dependent with respect to the time of symptom onset
    mutate(onset_known = as.logical(
      sapply(onset_time, function(t) rbernoulli(1, onset_known_prob[t]))
    ))

  reference_time <- ifelse(is.null(symp_to_rep_dist), "hosp_time", "onset_time")
  individual_cases$rep_time <- NA
  individual_cases$rep_time <- sapply(
    individual_cases[reference_time] %>% pull(),
    # rcat starts at 1, so need to shift
    function(t) t + extraDistr::rcat(1, delay_dist[t, ]) - 1
  )

  return(individual_cases)
}

#' Summarize information in simulated line list to various case counts by date.
summarize_linelist <- function(individual_cases, all_days) {
  onsets <- individual_cases %>%
    count(onset_time) %>%
    complete(onset_time = all_days, fill = list(n = 0)) %>%
    pull(n)
  infections_hospitalized <- individual_cases %>%
    count(infection_time) %>%
    complete(infection_time = all_days, fill = list(n = 0)) %>%
    pull(n)
  onsets_hospitalized <- individual_cases %>%
    count(onset_time) %>%
    complete(onset_time = all_days, fill = list(n = 0)) %>%
    pull(n)
  onsets_observed <- individual_cases %>%
    filter(onset_known) %>%
    count(onset_time) %>%
    complete(onset_time = all_days, fill = list(n = 0)) %>%
    pull(n)
  onsets_unobserved <- individual_cases %>%
    filter(!onset_known) %>%
    count(onset_time) %>%
    complete(onset_time = all_days, fill = list(n = 0)) %>%
    pull(n)
  hospitalizations <- individual_cases %>%
    count(hosp_time) %>%
    complete(hosp_time = all_days, fill = list(n = 0)) %>%
    pull(n)

  return(list(
    t = all_days,
    infections_hospitalized = infections_hospitalized[all_days],
    onsets_hospitalized = onsets_hospitalized[all_days],
    onsets_observed = onsets_observed[all_days],
    onsets_unobserved = onsets_unobserved[all_days],
    hospitalizations = hospitalizations[all_days]
  ))
}


#' Compute a mixture of discrete probability distributions.
#'
#' @param pmf_list A `list` of different probability mass functions.
#' @param weights A `vector` of same length as `pmf_list` with corresponding
#'   weights.
#'
#' @return A `vector` with the mixture PMF.
get_mixture_discrete <- function(pmf_list, weights) {
  l <- length(pmf_list[[1]])
  stopifnot(all(sapply(pmf_list, length) == l))
  stopifnot(length(weights) == length(pmf_list))
  weights <- weights / sum(weights)
  mix <- sapply(1:l, function(x) {
    sum(sapply(1:length(pmf_list), function(k) weights[k] * pmf_list[[k]][x]))
  })
  return(mix)
}

#' Compute a time series of distributions (one distribution per time step)
#' 
#' Discrete distributions are specified at certain changepoints. Between the
#' changepoints, the distribution transitions from the previous to the next
#' distribution using a mixture with weight parameter linearly or logistically
#' increasing from 0 to 1.
#'
#' @param changepoints A `vector` indicating different changepoints/
#' @param lengths As an alternative to `changepoints`, the length of each
#'   segment between two changepoints can be specified.
#' @param pmf_list A `list` with the discrete distributions at the different
#'   change points represented by vectors.
#' @param tribble_format As an alternative spceifying `changepoints` or
#'   `lengths` together with the `pmf_list`, a third option is to specify a
#'   tibble with a column for either `changepoints` or `lengths`, and a column
#'   with the `pmf_list`. This allows for a more readable specification, because
#'   the tibble can be defined conventiently using the `tribble` function, where
#'   changepoints and corresponding pmfs are directly next to each other.
#' @param shape Should the transition be "linear" or "logistic"
get_distribution_time_series <- function(
    changepoints = NULL, lengths = NULL,
    tribble_format = NULL, pmf_list, shape = "linear") {
  if (!is.null(tribble_format)) {
    if ("changepoints" %in% names(tribble_format)) {
      changepoints <- tribble_format$changepoints
    }
    if ("lengths" %in% names(tribble_format)) {
      lengths <- tribble_format$lengths
    }
    pmf_list <- tribble_format$pmf_list
  }
  if (is.null(changepoints)) {
    stopifnot(!is.null(lengths))
    if (length(lengths) == length(pmf_list)) {
      pmf_list[[length(lengths) + 1]] <- pmf_list[[length(lengths)]]
    }
    changepoints <- cumsum(c(1, lengths))
  }
  stopifnot(changepoints[1] == 1 & !is.unsorted(changepoints))
  start <- list(matrix(pmf_list[[1]], nrow = 1))
  transitions <- lapply(1:(length(changepoints) - 1), function(cp_index) {
    if (shape == "linear") {
      # weights for new dist grow from 0 to 1 according to a linear function
      mixing_weights <- seq(0, 1,
        length.out =
          changepoints[cp_index + 1] - changepoints[cp_index] + 1
      )[-1]
    } else if (shape == "logistic") {
      # weights for new dist grow from 0 to 1 according to a logistic function
      mixing_weights <- round(
        plogis(seq(-6, 6, length.out = changepoints[cp_index + 1] -
          changepoints[cp_index] + 1)), 2
      )[-1]
    } else {
      stop("No valid shape function ('linear' or 'logistic') provided.")
    }
    transition <- t(sapply(
      mixing_weights,
      function(w) {
        get_mixture_discrete(
          list(
            pmf_list[[cp_index]],
            pmf_list[[cp_index + 1]]
          ), c(1 - w, w)
        )
      }
    ))
    return(transition)
  })

  return(do.call(rbind, append(start, transitions)))
}

#' Get a piecewise stepwise time series.
#'
#' @param pieces Time points of the step function change points.
#' @param lengths Alternatively, the lengths of the segments between the change
#'   points.
#' @param values The different values on each segment.
#'
#' @return A `vector` with the piecewiese time series. The values are repeated
#'   during each segment and then jump to the next value at a changepoint
#'   according to a step function.
get_piecewise_time_series <- function(pieces, lengths = NULL, values) {
  if (is.null(lengths)) {
    lengths <- diff(pieces)
  }
  return(rep(values, times = lengths))
}

#' Create a time series with step function-like jumps between different scalar
#' values, and potentially seasonality added.
#'
#' @param changepoints A `vector` indicating different changepoints/
#' @param lengths As an alternative to `changepoints`, the length of each
#'   segment between two changepoints can be specified.
#' @param scalars A `vector` with the different values at the
#'   change points.
#' @param tribble_format As an alternative spceifying `changepoints` or
#'   `lengths` together with the `scalars`, a third option is to specify a
#'   tibble with a column for either `changepoints` or `lengths`, and a column
#'   with the `scalars`. This allows for a more readable specification, because
#'   the tibble can be defined conventiently using the `tribble` function, where
#'   changepoints and corresponding values are directly next to each other.
#' @param seasonality Should a seasonal variation / periodicity be added to the
#'   scalar time series.
get_scalar_time_series <- function(
    changepoints = NULL, lengths = NULL,
    tribble_format = NULL, scalars, seasonality = 0) {
  if (!is.null(tribble_format)) {
    if ("changepoints" %in% names(tribble_format)) {
      changepoints <- tribble_format$changepoints
    }
    if ("lengths" %in% names(tribble_format)) {
      lengths <- tribble_format$lengths
    }
    scalars <- tribble_format$scalars
  }
  if (is.null(changepoints)) {
    stopifnot(!is.null(lengths))
    if (length(lengths) == length(scalars)) {
      scalars <- c(scalars, scalars[length(lengths)])
    }
    changepoints <- cumsum(c(1, lengths))
  }
  stopifnot(changepoints[1] == 1 & !is.unsorted(changepoints))
  transitions <- sapply(1:(length(changepoints) - 1), function(cp_index) {
    return(seq(scalars[cp_index],
      scalars[cp_index + 1],
      length.out = changepoints[cp_index + 1] -
        changepoints[cp_index] + 1
    )[-1])
  })
  ts <- c(scalars[1], unlist(transitions))
  ts <- ts + rep_len(seasonality, length(ts))
  return(ts)
}


#' Get a time series of hazards (specifying a probability distribution) with
#' seasonality and trends.
#'
#' @param starting_hazard_logit A `vector` describing the baseline hazard on the
#'   logit scale. For example, this could be the hazard of reporting at
#'   different delays (starting with a delay of zero days).
#' @param changepoints A `vector` indicating different changepoints/
#' @param lengths As an alternative to `changepoints`, the length of each
#'   segment between two changepoints can be specified.
#' @param scalars A `vector` with different modifiers of the logit hazard
#'   according to a proportional hazards model.
#' @param tribble_format As an alternative spceifying `changepoints` or
#'   `lengths` together with the `scalars`, a third option is to specify a
#'   tibble with a column for either `changepoints` or `lengths`, and a column
#'   with the `scalars`. This allows for a more readable specification, because
#'   the tibble can be defined conventiently using the `tribble` function, where
#'   changepoints and corresponding values are directly next to each other.
#' @param ref_seasonality Seasonal variation by reference date (applied using a
#'   proportional hazards model).
#' @param rep_seasonality Seasonal variation by reporting date (applied using a
#'   proportional hazards model).
#'   
#' @return A matrix representing the hazard at each time step.
get_hazard_time_series <- function(
    starting_hazard_logit, changepoints = NULL, lengths = NULL,
    tribble_format = NULL, scalars, ref_seasonality = 0, rep_seasonality = 0) {
  get_deltas <- get_scalar_time_series(
    changepoints, lengths, tribble_format, scalars, ref_seasonality
  )
  hazard_logit_ts <- matrix(
    rep(starting_hazard_logit, length(get_deltas)),
    nrow = length(get_deltas),
    byrow = T
  )
  # add reference date trend and seasonality
  for (j in 1:ncol(hazard_logit_ts)) {
    hazard_logit_ts[, j] <- hazard_logit_ts[, j] + get_deltas
  }
  # add reporting date seasonality
  seasonality_len <- length(rep_seasonality)
  rep_seasonality_extd <- c(rep_seasonality, rep_seasonality)
  for (i in 1:nrow(hazard_logit_ts)) {
    seasonality_start <- 1 + ((i - 1) %% seasonality_len)
    seasonality_curr <- rep_seasonality_extd[
      seasonality_start:(seasonality_start + seasonality_len - 1)
    ]
    hazard_logit_ts[i, ] <- hazard_logit_ts[i, ] +
      rep_len(seasonality_curr, ncol(hazard_logit_ts))
  }
  # convert to unit scale
  hazard_ts <- plogis(hazard_logit_ts)
  return(hazard_ts)
}

#' Same as `get_hazard_time_series()`, but return probabilities instead of 
#' hazards.
#' 
#' @inheritParams get_hazard_time_series
#'
#' @return A matrix representing the probabilities at each time step.
get_p_time_series <- function(
    starting_hazard_logit, changepoints = NULL, lengths = NULL,
    tribble_format = NULL, scalars, ref_seasonality = 0, rep_seasonality = 0) {
  return(t(apply(get_hazard_time_series(
    starting_hazard_logit, changepoints, lengths,
    tribble_format, scalars, ref_seasonality, rep_seasonality
  ), 1, get_p_from_hazard)))
}

#' Add seasonal variation to a time series of probability distributions
#' 
#' The seasonality is added by converting the probability distributions to the
#' hazard scale, applying a proportional hazards model, and the converting back
#' to probabilities.
add_seasonality <- function(p_ts, rep_seasonality = 0) {
  # convert to hazard logit scale
  hazard_logit_ts <- qlogis(t(apply(p_ts, 1, get_hazard_from_p)))
  # add reporting date seasonality
  seasonality_len <- length(rep_seasonality)
  rep_seasonality_extd <- c(rep_seasonality, rep_seasonality)
  for (i in 1:nrow(hazard_logit_ts)) {
    seasonality_start <- 1 + ((i - 1) %% seasonality_len)
    seasonality_curr <- rep_seasonality_extd[
      seasonality_start:(seasonality_start + seasonality_len - 1)
    ]
    hazard_logit_ts[i, ] <- hazard_logit_ts[i, ] +
      rep_len(seasonality_curr, ncol(hazard_logit_ts))
  }
  # convert back to probabilities
  hazard_ts <- t(apply(plogis(hazard_logit_ts), 1, get_p_from_hazard))
  return(hazard_ts)
}

#' Helper function to conveniently plot a time series of scalar values.
plot_scalar_ts <- function(ts, x = NULL, xlab = "Time", ylab = "Value") {
  if (is.null(x)) x <- 1:length(ts)
  return(qplot(y = ts, x = x, geom = "line") +
    xlab(xlab) +
    ylab(ylab) +
    theme_bw())
}

#' Helper function to conveniently plot a time series of discrete probability
#' distributions.
plot_dist_ts <- function(ts_matrix, xlab = "Time", ylab = "Delay") {
  return(
    data.frame(ts_matrix) %>%
      setNames(1:ncol(ts_matrix)) %>%
      mutate(t = 1:n()) %>%
      pivot_longer(-t, names_to = "s", values_to = "p") %>%
      mutate(s = as.integer(s)) %>%
      ggplot(aes(t, s)) +
      geom_raster(aes(fill = p), interpolate = TRUE) +
      theme_bw() +
      coord_cartesian(expand = F) +
      xlab(xlab) +
      ylab(ylab)
  )
}

#' Helper function to conveniently plot quantiles of a time series of discrete
#' probability distributions.
plot_dist_ts_quantiles <- function(x) {
  df <- data.frame(
    mean = dist_get_mean(x),
    median = dist_get_quantile(x, p = 0.5),
    lower = dist_get_quantile(x, p = 0.05),
    upper = dist_get_quantile(x, p = 0.95)
  ) %>%
    mutate(t = 1:n())
  ggplot(df, aes(x = t)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3) +
    geom_line(aes(y = median)) +
    geom_line(aes(y = mean), linetype = "dashed") +
    ylab("x") +
    theme_bw()
}

#' Helper function to get the arithmetic mean of a matrix representing a
#' distribution over time.
dist_get_mean <- function(ts_matrix) {
  return(rowSums(ts_matrix * matrix(
    rep(0:(ncol(ts_matrix) - 1), nrow(ts_matrix)),
    ncol = ncol(ts_matrix), byrow = T
  )))
}

#' Helper function to get the median of a matrix representing a distribution
#' over time.
dist_get_quantile <- function(x, p = 0.5) {
  sapply(1:nrow(x), function(i) which(cumsum(x[i, ]) >= p)[1])
}
