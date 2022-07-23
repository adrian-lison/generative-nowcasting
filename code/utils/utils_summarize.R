
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
  fit_summary[["fit_datetime"]] <- fitted_model$metadata()$start_datetime
  fit_summary[["fit_seed"]] <- fitted_model$metadata()$seed
  fit_summary[["fit_duration"]] <- fitted_model$metadata()$time

  fit_summary[["start_date"]] <- start_date
  fit_summary[["now"]] <- now
  fit_summary[["T"]] <- stan_data_list[["T"]]
  fit_summary[["D"]] <- stan_data_list[["D"]]
  
  if (model_type == "impute") {
    fit_summary[["mu"]] <- fitted_model$summary("mu", c("median", "quantile2")) %>% 
    transmute(date = (1:n())+D, mu = median, .lower = q5, .upper = q95, .width = factor(0.95)) %>% 
    mutate(date = recode(date, !!!seq(start_date, now, by = 1)))
    
    fit_summary[["imputed"]] <- fitted_model %>%
      spread_draws(reported_unknown_imputed[date, d]) %>% 
      group_by(.draw, date) %>% 
      summarize(total = sum(reported_unknown_imputed), .groups = "drop") %>%
      group_by(date) %>% 
      median_qi(total, .width = c(0.5, 0.95)) %>%
      mutate(.width = as.factor(.width)) %>%
      select(-c(.point, .interval)) %>%
      mutate(date = date + stan_data_list[["D"]]) %>% # first D dates are not imputed
      mutate(date = recode(date, !!!seq(start_date, now, by = 1)))
    
    fit_summary[["imputed_posterior"]] <- get_imputed_posterior(fitted_model, stan_data_list[["D"]])
    
  } else {
    fit_summary[["L"]] <- stan_data_list[["L"]]
    fit_summary[["n_lambda_pre"]] <- stan_data_list[["n_lambda_pre"]]
    fit_summary[["max_gen"]] <- stan_data_list[["max_gen"]]
    
    fit_summary[["nowcast"]] <- summarize_nowcast(fitted_model, start_date, now)
    fit_summary[["delays"]] <- summarize_p(fitted_model, start_date, now)
    fit_summary[["mean_delay"]] <- summarize_mean_delay(fitted_model, start_date, now)
    fit_summary[["fraction_complete"]] <- summarize_alpha(fitted_model, start_date, now)
    
    if (model_type == "latent" | model_type == "renewal") {
      fit_summary[["latent"]] <- summarize_latent(fitted_model, start_date, now, fit_summary[["n_lambda_pre"]], fit_summary[["L"]])
      fit_summary[["median_incubation"]] <- which(cumsum(stan_data_list[["latent_delay_dist"]]) > 0.5)[1]
    }
    
    if (model_type == "latent") {
      fit_summary[["R"]] <- summarize_R_epiestim(fitted_model, start_date, now,
                                                 maxT = fit_summary[["T"]],
                                                 n_lambda_pre = fit_summary[["n_lambda_pre"]],
                                                 L = fit_summary[["L"]],
                                                 n_samples = 1000,
                                                 estimation_window = 3,
                                                 mean_serial_interval = 4.8,
                                                 std_serial_interval = 2.3,
                                                 mean_Re_prior = 1
      )
    }
    
    if (model_type == "renewal") {
      fit_summary[["R"]] <- summarize_R(fitted_model, start_date, now, fit_summary[["n_lambda_pre"]], fit_summary[["L"]], fit_summary[["max_gen"]])
      fit_summary[["median_generation"]] <- which(cumsum(stan_data_list[["generation_time_dist"]]) > 0.5)[1]
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
        gather_draws(`(alpha_logit_start)|(alpha_logit_sd)|(iota_log_ar_start)|(iota_log_ar_sd)|(R_level_start)|(R_trend_start)|(R_sd)|(ets_alpha)|(ets_beta)|(ets_phi)`, regex = TRUE) %>%
        median_qi(.width = 0.95) %>%
        select(-c(.width, .point, .interval))
    } else {
      fit_summary[["other"]] <- fitted_model %>%
        gather_draws(`(alpha_logit_start)|(alpha_logit_sd)|(lambda_log_level_start)|(lambda_log_trend_start)|(lambda_log_sd)|(ets_alpha)|(ets_beta)|(ets_phi)`, regex = TRUE) %>%
        median_qi(.width = 0.95) %>%
        select(-c(.width, .point, .interval))
    }
  }
  
  return(fit_summary)
}

## ----------------------------------------------------------------
##                  Individual summary functions                 -
## ----------------------------------------------------------------

summarize_nowcast <- function(fit, start_date, now) {
  if ("nowcast_unknown" %in% fit$metadata()$stan_variables) {
    regex_args <- quo(`(nowcast_known)|(nowcast_unknown)|(predicted_missing_rep)`[date])
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
    spread_draws(p_log[d, t]) %>%
    mutate(p = exp(p_log), d = d - 1) %>%
    select(-p_log) %>%
    median_qi(.width = 0.95) %>%
    select(-c(.width, .point, .interval)) %>%
    mutate(t = recode(t, !!!seq(start_date, now, by = 1))) %>%
    mutate(t_d = t + d - 1)
}

summarize_mean_delay <- function(fit, start_date, now) {
  fit %>%
    spread_draws(p_log[d, date]) %>%
    mutate(p = exp(p_log), d = d - 1) %>%
    select(-p_log) %>%
    group_by(.draw, date) %>%
    summarize(mean_delay = weighted.mean(d, w = p), .groups = "drop") %>%
    group_by(date) %>%
    median_qi(.width = 0.95) %>%
    select(-c(.width, .point, .interval)) %>%
    mutate(date = recode(date, !!!seq(start_date, now, by = 1)))
}

summarize_alpha <- function(fit, start_date, now) {
  fit %>%
    spread_draws(alpha_log[date]) %>%
    mutate(alpha = exp(alpha_log)) %>%
    select(-alpha_log) %>%
    median_qi(.width = 0.95) %>%
    select(-c(.width, .point, .interval)) %>%
    mutate(date = recode(date, !!!seq(start_date, now, by = 1)))
}

summarize_latent <- function(fit, start_date, now, n_lambda_pre, L) {
  fit %>%
    spread_draws(iota[date]) %>%
    median_qi(.width = 0.95) %>%
    select(-c(.width, .point, .interval)) %>%
    mutate(date = recode(date, !!!seq(as.Date(start_date) - n_lambda_pre - L, as.Date(now), by = 1)))
}

summarize_R <- function(fit, start_date, now, n_lambda_pre, L, max_gen) {
  fit %>%
    spread_draws(R[date]) %>%
    median_qi(.width = 0.95) %>%
    select(-c(.width, .point, .interval)) %>%
    mutate(date = recode(date, !!!seq(as.Date(start_date) - n_lambda_pre - L + max_gen, as.Date(now), by = 1)))
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

summarize_R_epiestim <- function(fit, start_date, now, maxT, n_lambda_pre, L, n_samples = 1000, estimation_window = 3, mean_serial_interval = 4.8, std_serial_interval = 2.3, mean_Re_prior = 1) {
  right_bound <- (maxT + L + n_lambda_pre) - (estimation_window - 1)
  t_start <- seq(L + n_lambda_pre, right_bound)
  t_end <- t_start + estimation_window - 1

  R_draws <- fit %>%
    spread_draws(iota[date]) %>%
    ungroup() %>%
    filter(.draw %in% sample(unique(.draw), n_samples)) %>%
    arrange(.draw, date) %>%
    group_by(.draw) %>%
    group_modify(~ .sample_R(.x$iota, mean_serial_interval, std_serial_interval, t_start, t_end, mean_Re_prior)) %>%
    transmute(date = recode(t, !!!seq(as.Date(start_date) - n_lambda_pre - L, as.Date(now), by = 1)), R)

  R_draws %<>%
    group_by(date) %>%
    median_qi(.width = 0.95) %>%
    select(-c(.width, .point, .interval))

  return(R_draws)
}

get_imputed_posterior <- function(fitted_model, D, ndraws = 100) {
  df <- fitted_model %>%
    spread_draws(reported_unknown_imputed[date, d], ndraws = ndraws) %>% 
    select(.draw, date, d, reported_unknown_imputed) %>% 
    mutate(date = date + D) %>% 
    arrange(.draw, date, d)
  
  .draw_unique <- unique(df$.draw)
  date_unique <- unique(df$date)
  d_unique <- unique(df$d)
  
  post_draws <- array(data = df$reported_unknown_imputed,
                      dim=c(length(.draw_unique), 
                            length(date_unique), 
                            length(d_unique)), 
                      dimnames=list(.draw_unique, date_unique, d_unique))
  
  return(post_draws)
}
