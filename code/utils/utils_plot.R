
###########################################################################
###########################################################################
###                                                                     ###
###                 PLOTTING FUNCTIONALITY FOR NOWCASTS                 ###
###                                                                     ###
###########################################################################
###########################################################################

plot_nowcast_all <- function(fit_summary, predLag = 1, time_window = 28, event1_axislabel = "Date", event2_axislabel = "Count") {
  now <- fit_summary[["now"]]
  mindate <- now - time_window # e.g. time_window=28 to show the last four weeks
  maxdate <- now - predLag

  max_cases <- fit_summary[["nowcast"]] %>%
    filter(.width == 0.95) %>%
    filter(date > mindate, date <= maxdate) %>%
    pull(nowcast_all.upper) %>%
    max()

  date_breaks <- seq.Date(as.Date(now) - 2, as.Date(mindate), by = -7)

  nowcast_plot <- fit_summary[["nowcast"]] %>%
    filter(.width == 0.95) %>%
    arrange(date) %>%
    mutate(
      is_weekend = lubridate::wday(date, label = TRUE) %in% c("Sat", "Sun"),
      nowcast_all_ma = rollapply(nowcast_all, 7, mean, align = "center", fill = NA)
    ) %>%
    filter(date > mindate, date <= maxdate) %>%
    ggplot(aes(x = date, y = nowcast_all)) +
    geom_tile(aes(
      y = 0,
      height = Inf,
      fill = is_weekend
    ), alpha = .3) +
    scale_fill_manual(values = c("NA", "grey"))

  if ("latent" %in% names(fit_summary)) {
    nowcast_plot <- nowcast_plot + geom_ribbon(data = fit_summary[["latent"]] %>% filter(date > mindate, date <= maxdate - fit_summary[["median_incubation"]]), aes(ymin = .lower, ymax = .upper, y = iota), fill = "#ffbf80", alpha = 0.3) +
      geom_line(data = fit_summary[["latent"]] %>% filter(date > mindate, date <= maxdate - fit_summary[["median_incubation"]]), aes(y = iota), color = "orange", linetype = "dashed")
  }

  (nowcast_plot +
    geom_ribbon(aes(ymin = nowcast_all.lower, ymax = nowcast_all.upper), fill = "lightblue", alpha = 0.6) + geom_line() + geom_line(aes(y = nowcast_all_ma), color = "#006600") +
    theme_bw() + ylab(event2_axislabel) + xlab(event1_axislabel) +
    scale_x_date(expand = c(0, 0), breaks = date_breaks, date_labels = "%b %d") +
    scale_y_continuous(expand = c(0, 0)) +
    theme(legend.position = "none") +
    coord_cartesian(ylim = c(0, max_cases * 1.05), xlim = c(mindate + 1, now)) +
    ggtitle(paste0("Nowcast from ", now, ": All cases (date estimated if missing)"))) %>%
    return()
}

plot_nowcast_observed <- function(fit_summary, predLag = 1, time_window = 28, reported_cases = NULL, event1_axislabel = "Date", event2_axislabel = "Count") {
  now <- fit_summary[["now"]]
  mindate <- now - time_window # e.g. time_window=28 to show the last four weeks
  maxdate <- now - predLag

  max_cases <- fit_summary[["nowcast"]] %>%
    filter(.width == 0.95) %>%
    filter(date > mindate, date <= maxdate) %>%
    pull(nowcast_known.upper) %>%
    max()

  date_breaks <- seq.Date(as.Date(now) - 2, as.Date(mindate), by = -7)

  plot <- fit_summary[["nowcast"]] %>%
    filter(.width == 0.95) %>%
    mutate(
      is_weekend = lubridate::wday(date, label = TRUE) %in% c("Sat", "Sun"),
      nowcast_known_ma = rollapply(nowcast_known, 7, mean, align = "center", fill = NA)
    ) %>%
    filter(date > mindate, date <= maxdate) %>%
    ggplot(aes(x = date, y = nowcast_known)) +
    geom_tile(aes(
      y = 0,
      height = Inf,
      fill = is_weekend
    ), alpha = .3) +
    scale_fill_manual(values = c("NA", "grey"))

  if (!is.null(reported_cases)) {
    plot <- plot + geom_bar(data = reported_cases %>%
      filter(event1_date > mindate, event1_date < now), aes(x = event1_date, y = n), stat = "identity", color = NA, fill = "grey")
  }

  plot <- plot +
    geom_ribbon(aes(ymin = nowcast_known.lower, ymax = nowcast_known.upper), fill = "lightblue", alpha = 0.6) + geom_line() + geom_line(aes(y = nowcast_known_ma), color = "red") +
    theme_bw() + ylab(event2_axislabel) + xlab(event1_axislabel) +
    scale_x_date(expand = c(0, 0), breaks = date_breaks, date_labels = "%b %d") +
    scale_y_continuous(expand = c(0, 0)) +
    theme(legend.position = "none") +
    coord_cartesian(ylim = c(0, max_cases * 1.05), xlim = c(mindate + 1, now)) +
    ggtitle(paste0("Nowcast from ", now, ": Only complete cases (no missing dates)"))

  plot %>% return()
}

plot_R <- function(fit_summary, predLag = 1, time_window = 28, event1_axislabel = "Date") {
  fit_summary[["R"]] %>%
    ggplot(aes(x = date, y = R)) +
    geom_ribbon(aes(ymin = .lower, ymax = .upper)) +
    geom_line()

  now <- fit_summary[["now"]]
  mindate <- now - time_window # e.g. time_window=28 to show the last four weeks
  maxdate <- now - predLag - fit_summary[["median_incubation"]]

  max_R <- fit_summary[["R"]] %>%
    filter(date > mindate, date <= maxdate) %>%
    pull(.upper) %>%
    max()

  date_breaks <- seq.Date(as.Date(now) - 2, as.Date(mindate), by = -7)

  plot <- fit_summary[["R"]] %>%
    mutate(is_weekend = lubridate::wday(date, label = TRUE) %in% c("Sat", "Sun")) %>%
    filter(date > mindate, date <= maxdate) %>%
    ggplot(aes(x = date, y = R)) +
    geom_tile(aes(
      y = 0,
      height = Inf,
      fill = is_weekend
    ), alpha = .3) +
    scale_fill_manual(values = c("NA", "grey")) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    geom_ribbon(aes(ymin = .lower, ymax = .upper), fill = "#79d279", alpha = 0.6) +
    geom_line() +
    theme_bw() +
    ylab("Effective reproduction number") +
    xlab(event1_axislabel) +
    scale_x_date(expand = c(0, 0), breaks = date_breaks, date_labels = "%b %d") +
    scale_y_continuous(expand = c(0, 0)) +
    theme(legend.position = "none") +
    coord_cartesian(ylim = c(0, max_R * 1.05), xlim = c(mindate + 1, now)) +
    ggtitle(paste0("Nowcast from ", now, ""))

  plot %>% return()
}

plot_delay <- function(fit_summary, event1_axislabel = "Date") {
  fit_summary[["mean_delay"]] %>%
    mutate(is_weekend = lubridate::wday(date, label = TRUE) %in% c("Sat", "Sun")) %>%
    ggplot(aes(x = date, y = mean_delay)) +
    geom_tile(aes(
      y = 3.2,
      height = Inf,
      fill = is_weekend
    ), alpha = .3) +
    scale_fill_manual(values = c("NA", "grey")) +
    geom_ribbon(aes(ymin = .lower, ymax = .upper), fill = "#aa80ff", alpha = 0.6) +
    geom_line() +
    theme_bw() +
    theme(legend.position = "none") +
    xlab(event1_axislabel) +
    ylab("Mean delay [days]") +
    scale_x_date(expand = c(0, 0), date_breaks = "2 weeks", date_labels = "%b %d") +
    ggtitle("Estimated mean reporting delay")
}
