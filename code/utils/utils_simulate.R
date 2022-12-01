#' Define simulation parameters
#'
#' @param initial_infections A vector as long as the maximum generation time with initial/seeding infections
#'
#' @param R A vector of the instantaneous effective reproduction number by date of infection
#'
#' @param generation_dist A vector of the generation time distribution. Must sum to one.
#'
#' @param incubation_dist A vector of the incubation time distribution. Must sum to one.
#'
#' @param hosp_prob A vector of the probability of hospitalization. by date of infection.
#'
#' @param symp_to_hosp_dist A vector of the time distribution from symptom onset to hospitalization. Must sum to one.
#'
#' @param onset_known_prob A vector of the probability of having the symptom onset recorded, by date of symptom onset.
#'
#' @param symp_to_rep_dist A vector of the time distribution from symptom onset to reporting. Must sum to one.
#'
#' @param hosp_to_rep_dist A vector of the time distribution from hospitalization to reporting. Must sum to one.
define_simulation <- function(initial_infections, R, generation_dist, incubation_dist, hosp_prob, symp_to_hosp_dist, onset_known_prob, symp_to_rep_dist, hosp_to_rep_dist) {
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
  infections <- simulate_infections(simulation_def$initial_infections, simulation_def$R, simulation_def$generation_dist)
  print(paste("Overall number of infections:", sum(infections[["I"]])))
  cases <- simulate_outcomes_cases(infections[["I"]], simulation_def$incubation_dist, simulation_def$hosp_prob, simulation_def$symp_to_hosp_dist)
  cases <- simulate_ascertainment_cases(cases, simulation_def$onset_known_prob, simulation_def$symp_to_rep_dist, simulation_def$hosp_to_rep_dist)
  process_summary <- summarize_linelist(cases, 1:length(infections[["I"]]))
  process_summary[["infections"]] <- infections[["I"]]
  process_summary[["infections_expected"]] <- infections[["I_expected"]]
  process_summary[["t_onset"]] <- 1:(length(infections[["I"]]) + ncol(simulation_def$incubation_dist) - 1)
  return(list(parameters = simulation_def, linelist = cases, process_summary = process_summary))
}

#' Simulate infections using a renewal process
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
#' @return A `list` with vectors for the expected infections and the realized
#' infections, from `1:length(R)` respectively.
#'
simulate_infections <- function(initial_infections, R, generation_dist) {
  maxGen <- length(generation_dist[1, ])
  assertthat::are_equal(length(initial_infections), maxGen)
  assertthat::are_equal(maxGen + length(R), dim(generation_dist)[1])

  I_expected <- rep(0, length(R))
  I <- c(initial_infections, rep(0, length(R)))
  for (t in 1:length(R)) {
    # renewal equation
    I_expected[t] <- R[t] * sum(
      I[(maxGen + t - 1):(maxGen + t - maxGen)] *
        generation_dist[cbind((maxGen + t - 1):(maxGen + t - maxGen), 1:maxGen)]
    )
    I[maxGen + t] <- rpois(1, I_expected[t])
  }
  return(list(I_expected = I_expected, I = I[(maxGen + 1):length(I)]))
}

simulate_outcomes_cases <- function(infections, incubation_dist, hosp_prob, symp_to_hosp_dist) {
  infections_hosp <- sapply(1:length(infections), function(t) rbinom(1, infections[t], hosp_prob[t])) # sample how many infections are hospitalized, assuming probability of hospitalization is time-dependent with respect to the time of infection
  onsets <- matrix(sapply(1:length(infections_hosp), function(t) rmultinom(1, infections_hosp[t], incubation_dist[t, ])), ncol = 1) # simulate incubation times
  onsets <- cbind(rep(1:length(infections_hosp), each = length(incubation_dist[1, ])), onsets) # add infection times
  onsets <- cbind(onsets[, 1], onsets[, 1] + rep(0:(length(incubation_dist[1, ]) - 1), times = length(infections_hosp)), onsets[, 2]) # compute onset times
  colnames(onsets) <- c("infection_time", "onset_time", "onset_count")
  onsets <- as.data.frame(onsets)

  stopifnot(length(hosp_prob) >= max(onsets$infection_time))
  stopifnot(dim(symp_to_hosp_dist)[1] >= max(onsets$onset_time))

  individual_cases <- onsets %>%
    uncount(onset_count) %>%
    # we assume that the distribution of the time from symptom onset to hospitalization is time-dependent with respect to the time of symptom onset
    mutate(hosp_time = sapply(onset_time, function(t) t + extraDistr::rcat(1, symp_to_hosp_dist[t, ])))

  return(individual_cases)
}

simulate_ascertainment_cases <- function(individual_cases, onset_known_prob, symp_to_rep_dist = NULL, hosp_to_rep_dist = NULL) {
  if (is.null(symp_to_rep_dist) & is.null(hosp_to_rep_dist)) stop("Either delay from symptom onset or from hospitalization must be provided.")
  stopifnot(length(onset_known_prob) >= max(individual_cases$onset_time))
  if (is.null(symp_to_rep_dist)) {
    delay_dist <- hosp_to_rep_dist
    stopifnot(dim(delay_dist)[1] >= max(individual_cases$hosp_time, na.rm = T))
  } else {
    delay_dist <- symp_to_rep_dist
    stopifnot(dim(delay_dist)[1] >= max(individual_cases$onset_time, na.rm = T))
  }
  individual_cases <- individual_cases %>%
    # we assume that the probability of the symptom onset being recorded is time-dependent with respect to the time of symptom onset
    mutate(onset_known = as.logical(sapply(onset_time, function(t) rbernoulli(1, onset_known_prob[t]))))

  reference_time <- ifelse(is.null(symp_to_rep_dist), "hosp_time", "onset_time")
  individual_cases$rep_time <- NA
  individual_cases$rep_time <- sapply(
    individual_cases[reference_time] %>% pull(),
    function(t) t + extraDistr::rcat(1, delay_dist[t, ]) - 1 # rcat starts at 1, so need to shift
  )

  return(individual_cases)
}

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

#' Compute a time series of distributions (one distribution per time step).
#' Discrete distributions are specified at certain changepoints. Between the
#' changepoints, the distribution transitions from the previous to the next
#' distribution using a mixture with weight parameter inearly increasing
#' from 0 to 1.
#'
get_distribution_time_series <- function(changepoints = NULL, lengths = NULL, tribble_format = NULL, pmf_list) {
  if (!is.null(tribble_format)) {
    if ("changepoints" %in% names(tribble_format)) changepoints <- tribble_format$changepoints
    if ("lengths" %in% names(tribble_format)) lengths <- tribble_format$lengths
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
    mixing_weights <- seq(0, 1, length.out = changepoints[cp_index + 1] - changepoints[cp_index] + 1)[-1]
    transition <- t(sapply(mixing_weights, function(w) get_mixture_discrete(list(pmf_list[[cp_index]], pmf_list[[cp_index + 1]]), c(1 - w, w))))
    return(transition)
  })

  return(do.call(rbind, append(start, transitions)))
}

get_piecewise_time_series <- function(pieces, lengths = NULL, values) {
  if (is.null(lengths)) {
    lengths <- diff(pieces)
  }
  return(rep(values, times = lengths))
}

get_scalar_time_series <- function(changepoints = NULL, lengths = NULL, tribble_format = NULL, scalars, seasonality = 0) {
  if (!is.null(tribble_format)) {
    if ("changepoints" %in% names(tribble_format)) changepoints <- tribble_format$changepoints
    if ("lengths" %in% names(tribble_format)) lengths <- tribble_format$lengths
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
    return(seq(scalars[cp_index], scalars[cp_index + 1], length.out = changepoints[cp_index + 1] - changepoints[cp_index] + 1)[-1])
  })
  ts <- c(scalars[1], unlist(transitions))
  ts <- ts + rep_len(seasonality, length(ts))
  return(ts)
}

get_hazard_time_series <- function(starting_hazard_logit, changepoints = NULL, lengths = NULL, tribble_format = NULL, scalars, ref_seasonality = 0, rep_seasonality = 0) {
  get_deltas <- get_scalar_time_series(changepoints, lengths, tribble_format, scalars, ref_seasonality)
  hazard_logit_ts <- matrix(rep(starting_hazard_logit, length(get_deltas)), nrow = length(get_deltas), byrow = T)
  # add reference date trend and seasonality
  for (j in 1:ncol(hazard_logit_ts)) {
    hazard_logit_ts[, j] <- hazard_logit_ts[, j] + get_deltas
  }
  # add reporting date seasonality
  seasonality_len <- length(rep_seasonality)
  rep_seasonality_extd <- c(rep_seasonality, rep_seasonality)
  for (i in 1:nrow(hazard_logit_ts)) {
    seasonality_start <- 1 + ((i - 1) %% seasonality_len)
    seasonality_curr <- rep_seasonality_extd[seasonality_start:(seasonality_start + seasonality_len - 1)]
    hazard_logit_ts[i, ] <- hazard_logit_ts[i, ] + rep_len(seasonality_curr, ncol(hazard_logit_ts))
  }
  # convert to unit scale
  hazard_ts <- plogis(hazard_logit_ts)
  return(hazard_ts)
}

add_seasonality <- function(p_ts, rep_seasonality = 0) {
  # convert to hazard logit scale
  hazard_logit_ts <- qlogis(t(apply(p_ts, 1, get_hazard_from_p)))
  # add reporting date seasonality
  seasonality_len <- length(rep_seasonality)
  rep_seasonality_extd <- c(rep_seasonality, rep_seasonality)
  for (i in 1:nrow(hazard_logit_ts)) {
    seasonality_start <- 1 + ((i - 1) %% seasonality_len)
    seasonality_curr <- rep_seasonality_extd[seasonality_start:(seasonality_start + seasonality_len - 1)]
    hazard_logit_ts[i, ] <- hazard_logit_ts[i, ] + rep_len(seasonality_curr, ncol(hazard_logit_ts))
  }
  # convert back to unit scale
  hazard_ts <- t(apply(plogis(hazard_logit_ts), 1, get_p_from_hazard))
  return(hazard_ts)
}

get_p_time_series <- function(starting_hazard_logit, changepoints = NULL, lengths = NULL, tribble_format = NULL, scalars, ref_seasonality = 0, rep_seasonality = 0) {
  return(t(apply(get_hazard_time_series(
    starting_hazard_logit, changepoints, lengths, tribble_format, scalars, ref_seasonality, rep_seasonality
  ), 1, get_p_from_hazard)))
}

plot_scalar_ts <- function(ts, xlab = "Time", ylab = "Value") {
  return(qplot(y = ts, x = 1:length(ts), geom = "line") +
    xlab(xlab) +
    ylab(ylab) +
    theme_bw())
}

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

sample_overall_delays <- function(t, samples = 10) {
  samples <- rcat(samples, prob = symp_to_hosp_dist[t, ]) + rcat(samples, prob = hosp_to_rep_dist[t, ])
  counts <- as.vector(table(cut(samples, breaks = 0:(maxSympToHosp + maxHospToRep + 1))))
  counts <- counts / sum(counts)
  return(counts)
}

dist_get_mean <- function(ts_matrix) {
  return(rowSums(ts_matrix * matrix(rep(0:(ncol(ts_matrix) - 1), nrow(ts_matrix)), ncol = ncol(ts_matrix), byrow = T)))
}
