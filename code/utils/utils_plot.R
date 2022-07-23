
###########################################################################
###########################################################################
###                                                                     ###
###                 PLOTTING FUNCTIONALITY FOR NOWCASTS                 ###
###                                                                     ###
###########################################################################
###########################################################################

plot_nowcast_observed <- function(fit_summary, cases_truth = NULL, reported_cases = NULL, predLag = 1, time_window = 28, breaks_resolution = NULL, show_weekends = T, event1_axislabel = "Date", event2_axislabel = "Count") {
  now <- fit_summary[["now"]]
  mindate <- now - time_window # e.g. time_window=28 to show the last four weeks
  maxdate <- now - predLag

  max_cases <- fit_summary[["nowcast"]] %>%
    filter(.width == 0.95) %>%
    filter(date > mindate, date <= maxdate) %>%
    pull(nowcast_known.upper) %>%
    max()

  plot <- fit_summary[["nowcast"]] %>%
    filter(.width == 0.95) %>%
    filter(date > mindate, date <= maxdate) %>%
    ggplot(aes(x = date, y = nowcast_known)) +
    geom_vline(xintercept = fit_summary$start_date, linetype = "dashed", color = "black") +
    geom_vline(xintercept = fit_summary$start_date + fit_summary$D, linetype = "dotted", color = "black") +
    {
      if (show_weekends) geom_tile(data = .getWeekends(mindate, maxdate), aes(y = 0, height = Inf, fill = is_weekend), alpha = .3)
    } +
    scale_fill_manual(values = c("NA", "grey")) +
    {
      if (!is.null(reported_cases)) geom_bar(data = reported_cases %>% filter(event1_date > mindate, event1_date < now), aes(x = event1_date, y = n), stat = "identity", color = NA, fill = "grey")
    } +
    {
      if (!is.null(cases_truth)) geom_line(data = cases_truth %>% filter(date > mindate, date < now), aes(y = cases), color = "red")
    } +
    geom_ribbon(aes(ymin = nowcast_known.lower, ymax = nowcast_known.upper), fill = "#b3d1ff", alpha = 0.6) +
    geom_line(color = "#0066ff") +
    theme_bw() +
    ylab(event2_axislabel) +
    xlab(event1_axislabel) +
    scale_x_date(expand = c(0, 0), breaks = .getDateBreaks(mindate, now - 2, breaks_resolution), date_labels = "%b %d") +
    scale_y_continuous(expand = c(0, 0)) +
    theme(legend.position = "none") +
    coord_cartesian(ylim = c(0, max_cases * 1.05), xlim = c(mindate + 1, now)) +
    ggtitle(paste0("Nowcast from ", now, ": Only complete cases (no missing dates)"))

  return(plot)
}

plot_nowcast_all <- function(fit_summary, cases_truth = NULL, infections_truth = NULL, show_infections = F, predLag = 1, time_window = 28, breaks_resolution = NULL, event1_axislabel = "Date", show_weekends = T, event2_axislabel = "Count") {
  now <- fit_summary[["now"]]
  mindate <- now - time_window # e.g. time_window=28 to show the last four weeks
  maxdate <- now - predLag

  max_cases <- fit_summary[["nowcast"]] %>%
    filter(.width == 0.95) %>%
    filter(date > mindate, date <= maxdate) %>%
    pull(nowcast_all.upper) %>%
    max()

  plot <- fit_summary[["nowcast"]] %>%
    filter(.width == 0.95) %>%
    arrange(date) %>%
    filter(date > mindate, date <= maxdate) %>%
    ggplot(aes(x = date, y = nowcast_all)) +
    geom_vline(xintercept = fit_summary$start_date, linetype = "dashed", color = "black") +
    geom_vline(xintercept = fit_summary$start_date + fit_summary$D, linetype = "dotted", color = "black") +
    {
      if (show_weekends) geom_tile(data = .getWeekends(mindate, maxdate), aes(y = 0, height = Inf, fill = is_weekend), alpha = .3)
    } +
    scale_fill_manual(values = c("NA", "grey")) +
    {
      if (show_infections & !is.null(infections_truth)) geom_line(data = infections_truth %>% filter(date > mindate, date <= maxdate - fit_summary[["median_incubation"]]), aes(y = infections), color = "red")
    } +
    {
      if (!is.null(cases_truth)) geom_line(data = cases_truth %>% filter(date > mindate, date <= maxdate), aes(y = cases), color = "red")
    } +
    {
      if (show_infections & "latent" %in% names(fit_summary)) geom_ribbon(data = fit_summary[["latent"]] %>% filter(date > mindate, date <= maxdate - fit_summary[["median_incubation"]]), aes(ymin = .lower, ymax = .upper, y = iota), fill = "#ffbf80", alpha = 0.3)
    } +
    {
      if (show_infections & "latent" %in% names(fit_summary)) geom_line(data = fit_summary[["latent"]] %>% filter(date > mindate, date <= maxdate - fit_summary[["median_incubation"]]), aes(y = iota), color = "orange", linetype = "dashed")
    } +
    geom_ribbon(aes(ymin = nowcast_all.lower, ymax = nowcast_all.upper), fill = "#b3d1ff", alpha = 0.6) +
    geom_line(color = "#0066ff") +
    theme_bw() +
    ylab(event2_axislabel) +
    xlab(event1_axislabel) +
    scale_x_date(expand = c(0, 0), breaks = .getDateBreaks(mindate, now - 2, breaks_resolution), date_labels = "%b %d") +
    scale_y_continuous(expand = c(0, 0)) +
    theme(legend.position = "none") +
    coord_cartesian(ylim = c(0, max_cases * 1.05), xlim = c(mindate + 1, now)) +
    ggtitle(paste0("Nowcast from ", now, ": All cases (date estimated if missing)"))

  return(plot)
}

plot_infections <- function(fit_summary, infections_truth = NULL, predLag = 1, time_window = 28, breaks_resolution = NULL, event1_axislabel = "Date", show_weekends = T, event2_axislabel = "Count") {
  now <- fit_summary[["now"]]
  mindate <- now - time_window # e.g. time_window=28 to show the last four weeks
  maxdate <- now - predLag

  max_infections <- fit_summary[["latent"]] %>%
    filter(date > mindate, date <= maxdate - fit_summary[["median_incubation"]]) %>%
    pull(.upper) %>%
    max()

  plot <- fit_summary[["latent"]] %>%
    arrange(date) %>%
    filter(date > mindate, date <= maxdate - fit_summary[["median_incubation"]]) %>%
    ggplot(aes(x = date, y = iota)) +
    geom_vline(xintercept = fit_summary$start_date, linetype = "dashed", color = "black") +
    geom_vline(xintercept = fit_summary$start_date + fit_summary$D, linetype = "dotted", color = "black") +
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
    coord_cartesian(ylim = c(0, max_infections * 1.05), xlim = c(mindate + 1, now)) +
    ggtitle(paste0("Nowcast from ", now, ": Infections"))

  return(plot)
}

plot_R <- function(fit_summary, R_truth = NULL, predLag = 1, time_window = 28, breaks_resolution = NULL, show_weekends = T, event1_axislabel = "Date") {
  now <- fit_summary[["now"]]
  mindate <- now - time_window # e.g. time_window=28 to show the last four weeks
  maxdate <- now - predLag

  max_R <- fit_summary[["R"]] %>%
    filter(date > mindate, date <= maxdate) %>%
    pull(.upper) %>%
    max()

  plot <- fit_summary[["R"]] %>%
    mutate(is_weekend = lubridate::wday(date, label = TRUE) %in% c("Sat", "Sun")) %>%
    filter(date > mindate, date <= maxdate) %>% # - fit_summary[["median_incubation"]]
    ggplot(aes(x = date, y = R)) +
    geom_vline(xintercept = fit_summary$start_date, linetype = "dashed", color = "black") +
    geom_vline(xintercept = fit_summary$start_date + fit_summary$D, linetype = "dotted", color = "black") +
    {
      if (show_weekends) geom_tile(data = .getWeekends(mindate, maxdate), aes(y = 0, height = Inf, fill = is_weekend), alpha = .3)
    } +
    scale_fill_manual(values = c("NA", "grey")) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    {
      if (!is.null(R_truth)) geom_line(data = R_truth, aes(x = date, y = R), color = "red")
    } +
    geom_ribbon(aes(ymin = .lower, ymax = .upper), fill = "#79d279", alpha = 0.6) +
    geom_line(color = "#006600") +
    theme_bw() +
    ylab("R") +
    xlab(event1_axislabel) +
    scale_x_date(expand = c(0, 0), breaks = .getDateBreaks(mindate, now - 2, breaks_resolution), date_labels = "%b %d") +
    scale_y_continuous(expand = c(0, 0), breaks = seq(0, 2.6, 0.2)) +
    theme(legend.position = "none") +
    coord_cartesian(ylim = c(0, max_R * 1.05), xlim = c(mindate + 1, now)) +
    ggtitle(paste0("Nowcast from ", now, ": Effective reproduction number"))

  return(plot)
}

plot_alpha <- function(fit_summary, alpha_truth = NULL, predLag = 1, time_window = 28, breaks_resolution = NULL, show_weekends = T, event1_axislabel = "Date") {
  now <- fit_summary[["now"]]
  mindate <- now - time_window # e.g. time_window=28 to show the last four weeks
  maxdate <- now - predLag

  max_alpha <- fit_summary[["fraction_complete"]] %>%
    filter(date > mindate, date <= maxdate) %>%
    pull(.upper) %>%
    max()

  plot <- fit_summary[["fraction_complete"]] %>%
    filter(date > mindate, date <= maxdate) %>%
    ggplot(aes(x = date, y = alpha)) +
    geom_vline(xintercept = fit_summary$start_date, linetype = "dashed", color = "black") +
    geom_vline(xintercept = fit_summary$start_date + fit_summary$D, linetype = "dotted", color = "black") +
    {
      if (show_weekends) geom_tile(data = .getWeekends(mindate, maxdate), aes(y = 0, height = Inf, fill = is_weekend), alpha = .3)
    } +
    scale_fill_manual(values = c("NA", "grey")) +
    geom_ribbon(aes(ymin = .lower, ymax = .upper), fill = "#ffccff", alpha = 0.6) +
    geom_line(color = "#ff66ff") +
    {
      if (!is.null(alpha_truth)) geom_line(data = alpha_truth, aes(x = date, y = alpha), color = "red")
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

plot_delay <- function(fit_summary, delay_truth = NULL, predLag = 1, time_window = 28, breaks_resolution = NULL, show_weekends = T, event1_axislabel = "Date") {
  now <- fit_summary[["now"]]
  mindate <- now - time_window # e.g. time_window=28 to show the last four weeks
  maxdate <- now - predLag

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
    geom_vline(xintercept = fit_summary$start_date, linetype = "dashed", color = "black") +
    geom_vline(xintercept = fit_summary$start_date + fit_summary$D, linetype = "dotted", color = "black") +
    {
      if (show_weekends) geom_tile(data = .getWeekends(mindate, maxdate), aes(y = 0, height = Inf, fill = is_weekend), alpha = .3)
    } +
    scale_fill_manual(values = c("NA", "grey")) +
    geom_ribbon(aes(ymin = .lower, ymax = .upper), fill = "#aa80ff", alpha = 0.6) +
    geom_line(color = "#9900cc") +
    {
      if (!is.null(delay_truth)) geom_line(data = delay_truth, aes(x = date, y = delay), color = "red")
    } +
    theme_bw() +
    ylab("Mean delay [days]") +
    xlab(event1_axislabel) +
    scale_x_date(expand = c(0, 0), breaks = .getDateBreaks(mindate, now - 2, breaks_resolution), date_labels = "%b %d") +
    scale_y_continuous(expand = c(0, 0)) +
    theme(legend.position = "none") +
    coord_cartesian(ylim = c(min_delay * 0.95, max_delay * 1.05), xlim = c(mindate + 1, now)) +
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
