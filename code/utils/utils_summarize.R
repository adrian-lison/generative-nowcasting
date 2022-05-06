
###########################################################################
###########################################################################
###                                                                     ###
###                 SUMMARY OF FITTED NOWCASTING MODELS                 ###
###                                                                     ###
###########################################################################
###########################################################################

## ----------------------------------------------------------------
##                    Overall summary function                   -
## ----------------------------------------------------------------

summarize_fit <- function(fitted_model,
                          model_type,
                          start_date,
                          now,
                          stan_data_list) {
  fit_summary <- list()

  fit_summary[["model"]] <- fitted_model$metadata()$model_name

  fit_summary[["start_date"]] <- start_date
  fit_summary[["now"]] <- now
  fit_summary[["T"]] <- stan_data_list[["T"]]
  fit_summary[["D"]] <- stan_data_list[["D"]]
  fit_summary[["L"]] <- stan_data_list[["L"]]

  fit_summary[["nowcast"]] <- summarize_nowcast(fitted_model, start_date, now)
  fit_summary[["delays"]] <- summarize_p(fitted_model, start_date, now)
  fit_summary[["mean_delay"]] <- summarize_mean_delay(fitted_model, start_date, now)
  fit_summary[["fraction_complete"]] <- summarize_alpha(fitted_model, start_date, now)

  if (model_type == "latent" | model_type == "renewal") {
    fit_summary[["latent"]] <- summarize_latent(fitted_model, start_date, now, fit_summary[["D"]], fit_summary[["L"]])
    fit_summary[["median_incubation"]] <- which(cumsum(rev(stan_data_list[["latent_delay_dist"]])) > 0.5)[1]
  }

  if (model_type == "latent") {
    fit_summary[["R"]] <- summarize_R_epiestim(fitted_model, start_date, now,
      maxT = fit_summary[["T"]],
      D = fit_summary[["D"]],
      L = fit_summary[["L"]],
      n_samples = 1000,
      estimation_window = 3,
      mean_serial_interval = 4.8,
      std_serial_interval = 2.3,
      mean_Re_prior = 1
    )
  }

  if (model_type == "renewal") {
    fit_summary[["R"]] <- summarize_R(fitted_model, start_date, now, fit_summary[["D"]], fit_summary[["L"]])
    fit_summary[["median_generation"]] <- which(cumsum(rev(stan_data_list[["generation_time_dist"]])) > 0.5)[1]
  }

  fit_summary[["nowcast_known_posterior"]] <- fitted_model %>%
    spread_draws(nowcast_known[date]) %>%
    select(-c(.chain, .iteration))

  if (model_type == "latent") {
    fit_summary[["other"]] <- fitted_model %>%
      gather_draws(alpha_logit_sd, iota_log_sd) %>%
      median_qi(.width = 0.95) %>%
      select(-c(.width, .point, .interval))
  } else if (model_type == "renewal") {
    fit_summary[["other"]] <- fitted_model %>%
      gather_draws(alpha_logit_sd, R_log_sd) %>%
      median_qi(.width = 0.95) %>%
      select(-c(.width, .point, .interval))
  } else {
    fit_summary[["other"]] <- fitted_model %>%
      gather_draws(alpha_logit_sd, lambda_log_sd) %>%
      median_qi(.width = 0.95) %>%
      select(-c(.width, .point, .interval))
  }

  fit_summary[["fit_datetime"]] <- fitted_model$metadata()$start_datetime
  fit_summary[["fit_seed"]] <- fitted_model$metadata()$seed
  fit_summary[["fit_duration"]] <- fitted_model$metadata()$time

  return(fit_summary)
}

## ----------------------------------------------------------------
##                  Individual summary functions                 -
## ----------------------------------------------------------------

summarize_nowcast <- function(fit, start_date, now) {
  if ("nowcast_unknown" %in% fit$metadata()$stan_variables) {
    regex_args <- quo(`(nowcast_known)|(nowcast_unknown)`[date])
  } else {
    regex_args <- quo(`(nowcast_known)`[date])
  }

  nowcast <- fit %>% spread_draws(!!regex_args, regex = TRUE)

  if ("nowcast_unknown" %in% fit$metadata()$stan_variables) {
    nowcast <- nowcast %>% mutate(nowcast_all = nowcast_known + nowcast_unknown)
  }

  nowcast %>%
    median_qi(.width = c(0.25, 0.5, 0.9, 0.95)) %>%
    mutate(.width = as.factor(.width)) %>%
    select(-c(.point, .interval)) %>%
    mutate(date = recode(date, !!!seq(start_date, now, by = 1))) %>%
    return()
}

summarize_p <- function(fit, start_date, now) {
  fit %>%
    spread_draws(p[t, d]) %>%
    median_qi(.width = 0.95) %>%
    select(-c(.width, .point, .interval)) %>%
    mutate(t = recode(t, !!!seq(start_date, now, by = 1))) %>%
    mutate(t_d = t + d - 1)
}

summarize_mean_delay <- function(fit, start_date, now) {
  fit %>%
    spread_draws(p[d, date]) %>%
    group_by(.draw, date) %>%
    summarize(mean_delay = weighted.mean(d, w = p), .groups = "drop") %>%
    group_by(date) %>%
    median_qi(.width = 0.95) %>%
    select(-c(.width, .point, .interval)) %>%
    mutate(date = recode(date, !!!seq(start_date, now, by = 1)))
}

summarize_alpha <- function(fit, start_date, now) {
  fit %>%
    spread_draws(alpha[date]) %>%
    median_qi(.width = 0.95) %>%
    select(-c(.width, .point, .interval)) %>%
    mutate(date = recode(date, !!!seq(start_date, now, by = 1)))
}

summarize_latent <- function(fit, start_date, now, D, L) {
  fit %>%
    spread_draws(iota[date]) %>%
    median_qi(.width = 0.95) %>%
    select(-c(.width, .point, .interval)) %>%
    mutate(date = recode(date, !!!seq(as.Date(start_date) - D - L, as.Date(now), by = 1)))
}

summarize_R <- function(fit, start_date, now, D, L) {
  fit %>%
    spread_draws(R[date]) %>%
    median_qi(.width = 0.95) %>%
    select(-c(.width, .point, .interval)) %>%
    mutate(date = recode(date, !!!seq(as.Date(start_date) - D - L, as.Date(now), by = 1)))
}

.sample_R <- function(incid, mean_serial_interval, std_serial_interval, t_start, t_end, mean_Re_prior = 1) {
  R_instantaneous <- estimate_R(
    incid = incid,
    method = "parametric_si",
    config = EpiEstim::make_config(
      list(
        mean_si = mean_serial_interval,
        std_si = std_serial_interval,
        t_start = t_start,
        t_end = t_end,
        mean_prior = mean_Re_prior
      )
    )
  )

  R_mean <- R_instantaneous$R$`Mean(R)`
  R_sd <- R_instantaneous$R$`Std(R)`
  # draw one sample from the posterior for R_t
  R_draw <- rgamma(length(R_mean), shape = (R_mean / R_sd)^2, rate = R_mean / (R_sd^2))
  return(data.frame(t = t_end, R = R_draw))
}

summarize_R_epiestim <- function(fit, start_date, now, maxT, D, L, n_samples = 1000, estimation_window = 3, mean_serial_interval = 4.8, std_serial_interval = 2.3, mean_Re_prior = 1) {
  right_bound <- (maxT + L + D) - (estimation_window - 1)
  t_start <- seq(L + D, right_bound)
  t_end <- t_start + estimation_window - 1

  R_draws <- fit %>%
    spread_draws(iota[date]) %>%
    ungroup() %>%
    filter(.draw %in% sample(unique(.draw), n_samples)) %>%
    arrange(.draw, date) %>%
    group_by(.draw) %>%
    group_modify(~ .sample_R(.x$iota, mean_serial_interval, std_serial_interval, t_start, t_end, mean_Re_prior)) %>%
    transmute(date = recode(t, !!!seq(as.Date(start_date) - D - L, as.Date(now), by = 1)), R)

  R_draws %<>%
    group_by(date) %>%
    median_qi(.width = 0.95) %>%
    select(-c(.width, .point, .interval))

  return(R_draws)
}
