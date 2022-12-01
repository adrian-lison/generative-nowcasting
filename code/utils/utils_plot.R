
###########################################################################
###########################################################################
###                                                                     ###
###                 PLOTTING FUNCTIONALITY FOR NOWCASTS                 ###
###                                                                     ###
###########################################################################
###########################################################################



plot_nowcast_observed <- function(fit_summary,
                                  truth = NULL,
                                  reported_cases = NULL,
                                  time_window = 28,
                                  breaks_resolution = NULL,
                                  show_weekends = T,
                                  show_data_horizons = T,
                                  y_lim = NULL,
                                  interval_width = 0.95,
                                  event1_axislabel = "Date",
                                  event2_axislabel = "Count") {
  now <- fit_summary[["now"]]
  mindate <- now - time_window # e.g. time_window=28 to show the last four weeks
  maxdate <- now

  max_cases <- fit_summary[["nowcast"]] %>%
    filter(.width == interval_width) %>%
    filter(date > mindate, date <= maxdate) %>%
    pull(nowcast_known.upper) %>%
    max()

  if (!is.null(truth) && ("delay" %in% names(truth))) {
    truth <- truth %>%
      filter(
        delay >= 0,
        delay <= fit_summary[["D"]]
      )
  }

  plot <- fit_summary[["nowcast"]] %>%
    filter(.width == interval_width) %>%
    filter(date > mindate, date <= maxdate) %>%
    ggplot(aes(x = date, y = nowcast_known)) +
    {
      if (show_data_horizons) geom_vline(xintercept = fit_summary$start_date, linetype = "dashed", color = "black")
    } +
    {
      if (show_data_horizons) geom_vline(xintercept = fit_summary$start_date + fit_summary$D, linetype = "dotted", color = "black")
    } +
    {
      if (show_weekends) geom_tile(data = .getWeekends(mindate, maxdate), aes(y = 0, height = Inf, fill = is_weekend), alpha = .3)
    } +
    scale_fill_manual(values = c("NA", "grey")) +
    {
      if (!is.null(reported_cases)) geom_col(data = reported_cases %>% filter(event1_date > mindate, event1_date < now), aes(x = event1_date, y = n), position = "identity", color = NA, fill = "grey")
    } +
    {
      if (!is.null(truth) && "cases_known" %in% names(truth)) geom_line(data = truth %>% filter(date > mindate, date < now), aes(y = cases_known), color = "red")
    } +
    geom_ribbon(aes(ymin = nowcast_known.lower, ymax = nowcast_known.upper), fill = "#b3d1ff", alpha = 0.6) +
    geom_line(color = "#0066ff") +
    theme_bw() +
    ylab(event2_axislabel) +
    xlab(event1_axislabel) +
    scale_x_date(expand = c(0, 0), breaks = .getDateBreaks(mindate, now - 2, breaks_resolution), date_labels = "%b %d") +
    scale_y_continuous(expand = c(0, 0)) +
    theme(legend.position = "none") +
    coord_cartesian(ylim = ifelse(rep(is.null(y_lim), 2), c(0, max_cases * 1.05), y_lim), xlim = c(mindate + 1, now)) +
    ggtitle(paste0("Nowcast from ", now, ": Only complete cases (no missing dates)"))

  return(plot)
}

plot_nowcast_all <- function(fit_summary,
                             truth = NULL,
                             infections_truth = NULL,
                             show_infections = F,
                             time_window = 28,
                             breaks_resolution = NULL,
                             show_weekends = T,
                             show_data_horizons = T,
                             y_lim = NULL,
                             interval_width = 0.95,
                             reported_cases = NULL, # unused argument, only for compatibility
                             event1_axislabel = "Date",
                             event2_axislabel = "Count") {
  now <- fit_summary[["now"]]
  mindate <- now - time_window # e.g. time_window=28 to show the last four weeks
  maxdate <- now

  max_cases <- fit_summary[["nowcast"]] %>%
    filter(.width == interval_width) %>%
    filter(date > mindate, date <= maxdate) %>%
    pull(nowcast_all.upper) %>%
    max()

  if (!is.null(truth) && ("delay" %in% names(truth))) {
    truth <- truth %>%
      filter(
        delay >= 0,
        delay <= fit_summary[["D"]]
      )
  }

  plot <- fit_summary[["nowcast"]] %>%
    filter(.width == interval_width) %>%
    arrange(date) %>%
    filter(date > mindate, date <= maxdate) %>%
    ggplot(aes(x = date, y = nowcast_all)) +
    {
      if (show_data_horizons) geom_vline(xintercept = fit_summary$start_date, linetype = "dashed", color = "black")
    } +
    {
      if (show_data_horizons) geom_vline(xintercept = fit_summary$start_date + fit_summary$D, linetype = "dotted", color = "black")
    } +
    {
      if (show_weekends) geom_tile(data = .getWeekends(mindate, maxdate), aes(y = 0, height = Inf, fill = is_weekend), alpha = .3)
    } +
    scale_fill_manual(values = c("NA", "grey")) +
    {
      if (show_infections && !is.null(infections_truth)) geom_line(data = infections_truth %>% filter(date > mindate, date <= maxdate - fit_summary[["median_incubation"]]), aes(y = infections), color = "red")
    } +
    {
      if (!is.null(truth) && "cases_all" %in% names(truth)) geom_line(data = truth %>% filter(date > mindate, date <= maxdate), aes(y = cases_all), color = "red")
    } +
    {
      if (show_infections && "latent" %in% names(fit_summary)) geom_ribbon(data = fit_summary[["latent"]] %>% filter(date > mindate, date <= maxdate - fit_summary[["median_incubation"]]), aes(ymin = .lower, ymax = .upper, y = iota), fill = "#ffbf80", alpha = 0.3)
    } +
    {
      if (show_infections && "latent" %in% names(fit_summary)) geom_line(data = fit_summary[["latent"]] %>% filter(date > mindate, date <= maxdate - fit_summary[["median_incubation"]]), aes(y = iota), color = "orange", linetype = "dashed")
    } +
    geom_ribbon(aes(ymin = nowcast_all.lower, ymax = nowcast_all.upper), fill = "#b3d1ff", alpha = 0.6) +
    geom_line(color = "#0066ff") +
    theme_bw() +
    ylab(event2_axislabel) +
    xlab(event1_axislabel) +
    scale_x_date(expand = c(0, 0), breaks = .getDateBreaks(mindate, now - 2, breaks_resolution), date_labels = "%b %d") +
    scale_y_continuous(expand = c(0, 0)) +
    theme(legend.position = "none") +
    coord_cartesian(ylim = ifelse(rep(is.null(y_lim), 2), c(0, max_cases * 1.05), y_lim), xlim = c(mindate + 1, now)) +
    ggtitle(paste0("Nowcast from ", now, ": All cases (date estimated if missing)"))

  return(plot)
}

plot_infections <- function(fit_summary,
                            infections_truth = NULL,
                            time_window = 28,
                            breaks_resolution = NULL,
                            event1_axislabel = "Date",
                            show_weekends = T,
                            show_data_horizons = T,
                            y_lim = NULL,
                            event2_axislabel = "Count") {
  now <- fit_summary[["now"]]
  mindate <- now - time_window # e.g. time_window=28 to show the last four weeks
  maxdate <- now

  max_infections <- fit_summary[["latent"]] %>%
    filter(date > mindate, date <= maxdate - fit_summary[["median_incubation"]]) %>%
    pull(.upper) %>%
    max()

  plot <- fit_summary[["latent"]] %>%
    arrange(date) %>%
    filter(date > mindate, date <= maxdate - fit_summary[["median_incubation"]]) %>%
    ggplot(aes(x = date, y = iota)) +
    {
      if (show_data_horizons) geom_vline(xintercept = fit_summary$start_date, linetype = "dashed", color = "black")
    } +
    {
      if (show_data_horizons) geom_vline(xintercept = fit_summary$start_date + fit_summary$D, linetype = "dotted", color = "black")
    } +
    {
      if (show_weekends) geom_tile(data = .getWeekends(mindate, maxdate), aes(y = 0, height = Inf, fill = is_weekend), alpha = .3)
    } +
    scale_fill_manual(values = c("NA", "grey")) +
    {
      if (!is.null(infections_truth)) geom_line(data = infections_truth %>% filter(date > mindate, date <= maxdate - fit_summary[["median_incubation"]]), aes(y = infections), color = "red")
    } +
    geom_ribbon(aes(ymin = .lower, ymax = .upper, y = iota), fill = "#ffbf80", alpha = 0.3) +
    geom_line(color = "orange", linetype = "dashed") +
    theme_bw() +
    ylab(event2_axislabel) +
    xlab(event1_axislabel) +
    scale_x_date(expand = c(0, 0), breaks = .getDateBreaks(mindate, now - 2, breaks_resolution), date_labels = "%b %d") +
    scale_y_continuous(expand = c(0, 0)) +
    theme(legend.position = "none") +
    coord_cartesian(ylim = ifelse(rep(is.null(y_lim), 2), c(0, max_infections * 1.05), y_lim), xlim = c(mindate + 1, now)) +
    ggtitle(paste0("Nowcast from ", now, ": Infections"))

  return(plot)
}

plot_R <- function(fit_summary,
                   truth = NULL,
                   time_window = 28,
                   breaks_resolution = NULL,
                   show_weekends = T,
                   show_data_horizons = T,
                   y_lim = NULL,
                   event1_axislabel = "Date of secondary infection",
                   reported_cases = NULL) {
  now <- fit_summary[["now"]]
  mindate <- now - time_window # e.g. time_window=28 to show the last four weeks
  maxdate <- now

  max_R <- fit_summary[["R"]] %>%
    filter(date > mindate, date <= maxdate) %>%
    pull(.upper) %>%
    max()

  plot <- fit_summary[["R"]] %>%
    mutate(is_weekend = lubridate::wday(date, label = TRUE) %in% c("Sat", "Sun")) %>%
    filter(date > mindate, date <= maxdate) %>% # - fit_summary[["median_incubation"]]
    ggplot(aes(x = date, y = R)) +
    {
      if (show_data_horizons) geom_vline(xintercept = fit_summary$start_date, linetype = "dashed", color = "black")
    } +
    {
      if (show_data_horizons) geom_vline(xintercept = fit_summary$start_date + fit_summary$D, linetype = "dotted", color = "black")
    } +
    {
      if (show_weekends) geom_tile(data = .getWeekends(mindate, maxdate), aes(y = 0, height = Inf, fill = is_weekend), alpha = .3)
    } +
    scale_fill_manual(values = c("NA", "grey")) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    {
      if (!is.null(truth) && ("R" %in% names(truth))) geom_line(data = truth, aes(x = date, y = R), color = "red")
    } +
    geom_ribbon(aes(ymin = .lower, ymax = .upper), fill = "#79d279", alpha = 0.6) +
    geom_line(color = "#006600") +
    theme_bw() +
    ylab("R") +
    xlab(event1_axislabel) +
    scale_x_date(expand = c(0, 0), breaks = .getDateBreaks(mindate, now - 2, breaks_resolution), date_labels = "%b %d") +
    scale_y_continuous(expand = c(0, 0), breaks = seq(0, 8, 0.2)) +
    theme(legend.position = "none") +
    coord_cartesian(ylim = ifelse(rep(is.null(y_lim), 2), c(0, max_R * 1.05), y_lim), xlim = c(mindate + 1, now)) +
    ggtitle(paste0("Nowcast from ", now, ": Effective reproduction number"))

  return(plot)
}

plot_alpha <- function(fit_summary,
                       truth = NULL,
                       time_window = 28,
                       breaks_resolution = NULL,
                       show_weekends = T,
                       show_data_horizons = T,
                       event1_axislabel = "Date",
                       reported_cases = NULL) {
  now <- fit_summary[["now"]]
  mindate <- now - time_window # e.g. time_window=28 to show the last four weeks
  maxdate <- now

  print(fit_summary[["fraction_complete"]])

  max_alpha <- fit_summary[["fraction_complete"]] %>%
    filter(date > mindate, date <= maxdate) %>%
    pull(.upper) %>%
    max()

  plot <- fit_summary[["fraction_complete"]] %>%
    filter(date > mindate, date <= maxdate) %>%
    ggplot(aes(x = date, y = alpha)) +
    {
      if (show_data_horizons) geom_vline(xintercept = fit_summary$start_date, linetype = "dashed", color = "black")
    } +
    {
      if (show_data_horizons) geom_vline(xintercept = fit_summary$start_date + fit_summary$D, linetype = "dotted", color = "black")
    } +
    {
      if (show_weekends) geom_tile(data = .getWeekends(mindate, maxdate), aes(y = 0, height = Inf, fill = is_weekend), alpha = .3)
    } +
    scale_fill_manual(values = c("NA", "grey")) +
    geom_ribbon(aes(ymin = .lower, ymax = .upper), fill = "#ffccff", alpha = 0.6) +
    geom_line(color = "#ff66ff") +
    {
      if (!is.null(truth) && ("alpha" %in% names(truth))) geom_line(data = truth, aes(x = date, y = alpha), color = "red")
    } +
    theme_bw() +
    ylab("Known fraction") +
    xlab(event1_axislabel) +
    scale_x_date(expand = c(0, 0), breaks = .getDateBreaks(mindate, now - 2, breaks_resolution), date_labels = "%b %d") +
    scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
    theme(legend.position = "none") +
    coord_cartesian(ylim = c(0, max_alpha * 1.05), xlim = c(mindate + 1, now)) +
    ggtitle(paste0("Nowcast from ", now, ": Fraction of known onset dates"))

  return(plot)
}

plot_delay <- function(fit_summary,
                       truth = NULL,
                       time_window = 28,
                       breaks_resolution = NULL,
                       show_weekends = T,
                       show_data_horizons = T,
                       reported_cases = NULL, # unused argument, only for compatibility
                       event1_axislabel = "Date") {
  now <- fit_summary[["now"]]
  mindate <- now - time_window # e.g. time_window=28 to show the last four weeks
  maxdate <- now

  max_delay <- fit_summary[["mean_delay"]] %>%
    filter(date > mindate, date <= maxdate) %>%
    pull(.upper) %>%
    max()

  min_delay <- fit_summary[["mean_delay"]] %>%
    filter(date > mindate, date <= maxdate) %>%
    pull(.lower) %>%
    min()

  plot <- fit_summary[["mean_delay"]] %>%
    filter(date > mindate, date <= maxdate) %>%
    ggplot(aes(x = date, y = mean_delay)) +
    {
      if (show_data_horizons) geom_vline(xintercept = fit_summary$start_date, linetype = "dashed", color = "black")
    } +
    {
      if (show_data_horizons) geom_vline(xintercept = fit_summary$start_date + fit_summary$D, linetype = "dotted", color = "black")
    } +
    {
      if (show_weekends) geom_tile(data = .getWeekends(mindate, maxdate), aes(y = 0, height = Inf, fill = is_weekend), alpha = .3)
    } +
    scale_fill_manual(values = c("NA", "grey")) +
    geom_ribbon(aes(ymin = .lower, ymax = .upper), fill = "#aa80ff", alpha = 0.6) +
    geom_line(color = "#9900cc") +
    {
      if (!is.null(truth)) geom_line(data = truth, aes(x = date, y = mean_delay), color = "red")
    } +
    theme_bw() +
    ylab("Mean delay [days]") +
    xlab(event1_axislabel) +
    scale_x_date(expand = c(0, 0), breaks = .getDateBreaks(mindate, now - 2, breaks_resolution), date_labels = "%b %d") +
    scale_y_continuous(expand = c(0, 0)) +
    theme(legend.position = "none") +
    coord_cartesian(ylim = c(min_delay * 0.95, max_delay * 1.05), xlim = c(mindate + 1, maxdate)) +
    ggtitle(paste0("Nowcast from ", now, ": Estimated mean reporting delay"))

  return(plot)
}

.getWeekends <- function(from, to) {
  data.frame(date = seq.Date(as.Date(from), as.Date(to), by = 1)) %>%
    mutate(is_weekend = lubridate::wday(date, label = TRUE) %in% c("Sat", "Sun")) %>%
    return()
}

.getDateBreaks <- function(from, to, breaks_resolution = NULL) {
  if (is.null(breaks_resolution)) {
    date_breaks <- waiver()
  } else {
    date_breaks <- seq.Date(as.Date(to), as.Date(from), by = -breaks_resolution)
  }
  return(date_breaks)
}

#' Helper function, provides an interface for validation result nowcasts to
#' be plotted using the other plotting functions
plot_validate <- function(nowcasts, maxDelay, delay_plot, plot_nowcast_date, plot_function, time_window, with_naive = F, ...) {
  if (missing(delay_plot)) delay_plot <- unique(nowcasts$delay)
  if (missing(plot_nowcast_date)) plot_nowcast_date <- unique(nowcasts$nowcast_date)

  plot_nowcast_date <- as.Date(plot_nowcast_date)

  if (with_naive) {
    rep_cases <- nowcasts %>%
      select(nowcast_date, event1_date = date, n = nowcast_known_naive) %>%
      unique() %>%
      complete(
        event1_date = seq.Date(min(event1_date, na.rm = T),
          max(nowcast_date, na.rm = T),
          by = "1 day"
        ),
        nowcast_date = seq.Date(min(nowcast_date, na.rm = T),
          max(nowcast_date, na.rm = T),
          by = "1 day"
        )
      ) %>%
      filter(nowcast_date >= event1_date) %>%
      arrange(event1_date, nowcast_date) %>%
      group_by(event1_date) %>%
      fill(n, .direction = "downup") %>%
      mutate(n = na.fill(n, 0)) %>%
      ungroup() %>%
      filter(nowcast_date <= max(plot_nowcast_date, na.rm = T), !is.na(n))
  } else {
    rep_cases <- NULL
  }

  if ("R" %in% names(nowcasts)) {
    input_R <- nowcasts %>%
      mutate(.width = 0.95) %>%
      filter(delay %in% delay_plot, nowcast_date %in% plot_nowcast_date) %>%
      select(date, delay, nowcast_date, R, .upper = R.upper, .lower = R.lower)
  } else {
    input_R <- NULL
  }

  if ("mean_delay" %in% names(nowcasts)) {
    input_mean_delay <- nowcasts %>%
      mutate(.width = 0.95) %>%
      filter(delay %in% delay_plot, nowcast_date %in% plot_nowcast_date) %>%
      select(date, delay, nowcast_date, mean_delay, .upper = mean_delay.upper, .lower = mean_delay.lower)
  } else {
    input_mean_delay <- NULL
  }

  if ("alpha" %in% names(nowcasts)) {
    input_alpha <- nowcasts %>%
      mutate(.width = 0.95) %>%
      filter(delay %in% delay_plot, nowcast_date %in% plot_nowcast_date) %>%
      select(date, delay, nowcast_date, alpha, .upper = alpha.upper, .lower = alpha.lower)
  } else {
    input_alpha <- NULL
  }

  if (missing(time_window)) time_window <- max(nowcasts$nowcast_date, na.rm = T) - min(nowcasts$date, na.rm = T)

  return(plot_function(
    fit_summary = list(
      now = max(nowcasts$nowcast_date, na.rm = T),
      nowcast = nowcasts %>% mutate(.width = 0.95) %>%
        filter(delay %in% delay_plot, nowcast_date %in% plot_nowcast_date),
      R = input_R,
      mean_delay = input_mean_delay,
      fraction_complete = input_alpha,
      D = maxDelay
    ),
    time_window = time_window,
    truth = nowcasts %>%
      select(
        date,
        delay,
        any_of(c(
          cases_known = "nowcast_known_true",
          cases_unknown = "nowcast_unknown_true",
          cases_all = "nowcast_all_true",
          R = "R_true",
          mean_delay = "mean_delay_true",
          alpha = "alpha_true"
        ))
      ),
    reported_cases = rep_cases,
    ...
  ))
}
