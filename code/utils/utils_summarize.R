
# Summary of fitted nowcasting models ----

draws_R_renewal <- function(p_nowcast, start_date, now, ndraws, samples_per_draw, standata, inits) {
  renewal_model <- cmdstan_model(here::here("code", "models", "renewal.stan"))

  Rdraws <- lapply(1:ndraws, function(i) {
    onsets <- p_nowcast %>%
      filter(.draw == i) %>%
      pull(nowcast_all)
    sample_R_renewal(onsets, standata, renewal_model, inits, samples_per_draw)
  })

  Rdraws <- bind_rows(Rdraws) %>%
    mutate(date = recode(date, !!!seq(as.Date(start_date) - standata$n_lambda_pre - standata$L + standata$max_gen, as.Date(now), by = 1)))

  return(Rdraws)
}

sample_R_renewal <- function(onsets, standata, renewal_model, inits, samples_per_draw) {
  standata$T <- length(onsets)
  standata$onsets <- onsets

  Rfit <- renewal_model$sample(
    data = standata,
    iter_warmup = 500,
    iter_sampling = as.integer(ceiling(samples_per_draw / 4)),
    adapt_delta = 0.95,
    step_size = 0.01,
    max_treedepth = 15,
    chains = 4,
    parallel_chains = 1,
    seed = 42,
    refresh = 200,
    init = inits,
    show_messages = F
  )

  Rdraws <- Rfit %>%
    spread_draws(R[date])

  return(Rdraws)
}

## Overall summary function ----

summarize_fit <- function(fitted_model,
                          output_def,
                          model_type,
                          start_date,
                          now,
                          stan_data_list) {
  sampling_summary <- list(
    model = fitted_model$metadata()$model_name,
    fit_datetime = fitted_model$metadata()$start_datetime,
    fit_seed = fitted_model$metadata()$seed,
    fit_duration = fitted_model$metadata()$time
  )

  model_summary <- summarize_model(start_date, now, stan_data_list)

  if (model_type %in% c("impute", "impute_parametric")) {
    fit_summary <- summarize_impute(fitted_model, output_def, start_date, now, stan_data_list)
  } else {
    fit_summary <- list()

    fit_summary[["nowcast"]] <- summarize_nowcast(fitted_model, start_date, now)

    fit_summary[["mean_delay"]] <- summarize_mean_delay(fitted_model, start_date, now)

    if (!(model_type %in% c("nowcast_imputed", "nowcast_imputed_renewal"))) {
      fit_summary[["fraction_complete"]] <- summarize_alpha(fitted_model, start_date, now)
    }

    if (model_type %in% c("latent", "renewal", "nowcast_imputed_renewal")) {
      fit_summary[["latent"]] <- summarize_latent(fitted_model, start_date, now, stan_data_list$n_lambda_pre, stan_data_list$L)
      fit_summary[["median_incubation"]] <- which(cumsum(stan_data_list$latent_delay_dist) > 0.5)[1]
    }

    posterior_nowcast_temp <- posterior_nowcast(fitted_model, start_date, now)

    if (("posterior_nowcast" %in% output_def)) {
      fit_summary[["posterior_nowcast"]] <- posterior_nowcast_temp
    }

    {
      if (model_type %in% c("base", "nowcast_imputed")) {

        # epiEstim
        R_draws_epiestim <- draws_R_epiestim(fitted_model, start_date, now,
          maxT = stan_data_list$T,
          L = stan_data_list$L,
          latent_delay_dist = stan_data_list$latent_delay_dist,
          n_samples = 1000,
          estimation_window = 7,
          mean_serial_interval = 4.8,
          std_serial_interval = 2.3,
          mean_Re_prior = 1,
          sample_incubation = FALSE
        )

        R_draws_renewal <- draws_R_renewal(posterior_nowcast_temp, start_date, now,
          ndraws = 40,
          samples_per_draw = 100,
          standata = stan_data_list,
          inits = default_inits(stan_data_list, "renewal")
        )

        fit_summary[["R_epiestim"]] <- R_draws_epiestim %>%
          group_by(date) %>%
          median_qi(.width = 0.95) %>%
          select(-c(.width, .point, .interval))

        fit_summary[["R_renewal"]] <- R_draws_renewal %>%
          group_by(date) %>%
          median_qi(.width = 0.95) %>%
          select(-c(.width, .point, .interval))

        if ("posterior_R" %in% output_def) {
          fit_summary[["posterior_R_epiestim"]] <- R_draws_epiestim
          fit_summary[["posterior_R_renewal"]] <- R_draws_renewal
        }
      }
    }

    if (model_type %in% c("renewal", "nowcast_imputed_renewal")) {
      fit_summary[["R_renewal"]] <- summarize_R(fitted_model, start_date, now, stan_data_list$n_lambda_pre, stan_data_list$L, stan_data_list$max_gen)
      fit_summary[["median_generation"]] <- which(cumsum(stan_data_list$generation_time_dist) > 0.5)[1]

      if ("posterior_R" %in% output_def) {
        fit_summary[["posterior_R_renewal"]] <- posterior_R(fitted_model, start_date, now, stan_data_list$n_lambda_pre, stan_data_list$L, stan_data_list$max_gen)
      }
    }

    if ("delays" %in% output_def) {
      fit_summary[["delays"]] <- summarize_p(fitted_model, start_date, now)
    }

    fit_summary[["other"]] <- fitted_model %>%
      gather_draws(`(alpha_logit_start)|(alpha_logit_sd)|(lambda_log_level_start)|(lambda_log_trend_start)|(lambda_log_sd)|(iota_log_ar_start)|(iota_log_ar_sd)|(iota_log_sd)|(R_level_start)|(R_trend_start)|(R_sd)|(ets_alpha)|(ets_beta)|(ets_phi)|(xi_negbinom)`, regex = TRUE) %>%
      median_qi(.width = 0.95) %>%
      select(-c(.width, .point, .interval))
  }

  return(c(sampling_summary, model_summary, fit_summary))
}

## Model summary functions ----

summarize_model <- function(start_date, now, stan_data_list) {
  model_summary <- list(
    start_date = start_date,
    now = now
  )
  for (v in c("T", "D", "L", "n_lambda_pre", "max_gen")) {
    if (v %in% names(stan_data_list)) {
      model_summary[[v]] <- stan_data_list[[v]]
    }
  }
  return(model_summary)
}

summarize_impute <- function(fitted_model,
                             output_def,
                             start_date,
                             now,
                             stan_data_list) {
  impute_summary <- list()

  if ("mu" %in% fitted_model$metadata()$stan_variables) {
    impute_summary[["mu"]] <- fitted_model$summary("mu", c("median", "quantile2")) %>%
      transmute(
        date = (1:n()) + stan_data_list[["D"]],
        mu = median,
        .lower = q5,
        .upper = q95,
        .width = factor(0.95)
      ) %>%
      mutate(date = recode(date, !!!seq(start_date, now, by = 1)))
  }

  if ("p" %in% fitted_model$metadata()$stan_variables) {
    impute_summary[["mean_delay_backward"]] <- summarize_mean_delay(fitted_model, start_date + stan_data_list$D, now)
  }

  impute_summary[["imputed"]] <- fitted_model %>%
    spread_draws(reported_unknown_imputed[date, d]) %>%
    group_by(.draw, date) %>%
    summarize(total = sum(reported_unknown_imputed), .groups = "drop") %>%
    group_by(date) %>%
    median_qi(total, .width = c(0.5, 0.95)) %>%
    mutate(.width = as.factor(.width)) %>%
    select(-c(.point, .interval)) %>%
    mutate(date = date + stan_data_list[["D"]]) %>% # first D dates are not imputed
    mutate(date = recode(date, !!!seq(start_date, now, by = 1)))

  impute_summary[["imputed_posterior"]] <- get_imputed_posterior(
    fitted_model,
    stan_data_list[["D"]],
    ndraws = 1
  )

  return(impute_summary)
}

## Individual summary functions ----

summarize_nowcast <- function(fit, start_date, now) {
  regex_args <- quo(`(nowcast_known)|(nowcast_unknown)|(nowcast_all)|(predicted_missing_rep)`[date])
  nowcast <- fit %>% spread_draws(!!regex_args, regex = TRUE)

  if (all(c("nowcast_known", "nowcast_unknown") %in% names(nowcast)) &
    !("nowcast_all" %in% names(nowcast))) {
    nowcast <- nowcast %>% mutate(nowcast_all = nowcast_known + nowcast_unknown)
  }

  nowcast %>%
    median_qi(.width = c(0.5, 0.95), .simple_names = F) %>%
    mutate(.width = as.factor(.width)) %>%
    select(-c(.point, .interval)) %>%
    mutate(date = recode(date, !!!seq(start_date, now, by = 1))) %>%
    return()
}

posterior_nowcast <- function(fit, start_date, now) {
  regex_args <- quo(`(nowcast_known)|(nowcast_unknown)|(nowcast_all)|(predicted_missing_rep)`[date])
  nowcast <- fit %>%
    spread_draws(!!regex_args, regex = TRUE) %>%
    select(-c(.chain, .iteration)) %>%
    mutate(date = recode(date, !!!seq(start_date, now, by = 1)))

  if (all(c("nowcast_known", "nowcast_unknown") %in% names(nowcast)) &
    !("nowcast_all" %in% names(nowcast))) {
    nowcast <- nowcast %>% mutate(nowcast_all = nowcast_known + nowcast_unknown)
  }

  return(nowcast)
}

summarize_p <- function(fit, start_date, now, by_reference = TRUE) {
  if ("p" %in% fitted_model$metadata()$stan_variables) {
    p_summary <- fit %>%
      spread_draws(p[d, t]) %>%
      mutate(d = d - 1) %>%
      median_qi(.width = 0.95) %>%
      select(-c(.width, .point, .interval)) %>%
      mutate(t = recode(t, !!!seq(start_date, now, by = 1)))
  } else if ("p_log" %in% fitted_model$metadata()$stan_variables) {
    p_summary <- fit %>%
      spread_draws(p_log[d, t]) %>%
      mutate(p = exp(p_log), d = d - 1) %>%
      select(-p_log) %>%
      median_qi(.width = 0.95) %>%
      select(-c(.width, .point, .interval)) %>%
      mutate(t = recode(t, !!!seq(start_date, now, by = 1)))
  } else {
    stop("No variable p or p_log found in fitted model.")
  }

  if (by_reference) {
    p_summary <- p_summary %>%
      mutate(t_d = t + (d - 1))
  } else {
    # by reporting date
    p_summary <- p_summary %>%
      mutate(t_d = t - (d - 1))
  }
  return(p_summary)
}

summarize_mean_delay <- function(fit, start_date, now) {
  if ("p" %in% fit$metadata()$stan_variables) {
    fit %>%
      spread_draws(p[d, date]) %>%
      mutate(d = d - 1) %>%
      group_by(.draw, date) %>%
      summarize(mean_delay = weighted.mean(d, w = p), .groups = "drop") %>%
      group_by(date) %>%
      median_qi(.width = 0.95) %>%
      select(-c(.width, .point, .interval)) %>%
      mutate(date = recode(date, !!!seq(start_date, now, by = 1))) %>%
      return()
  } else if ("p_log" %in% fit$metadata()$stan_variables) {
    fit %>%
      spread_draws(p_log[d, date]) %>%
      mutate(p = exp(p_log), d = d - 1) %>%
      select(-p_log) %>%
      group_by(.draw, date) %>%
      summarize(mean_delay = weighted.mean(d, w = p), .groups = "drop") %>%
      group_by(date) %>%
      median_qi(.width = 0.95) %>%
      select(-c(.width, .point, .interval)) %>%
      mutate(date = recode(date, !!!seq(start_date, now, by = 1))) %>%
      return()
  } else {
    stop("No variable p or p_log found in fitted model.")
  }
}

summarize_alpha <- function(fit, start_date, now) {
  fit %>%
    spread_draws(alpha_log[date]) %>%
    mutate(alpha = exp(alpha_log)) %>%
    select(-alpha_log) %>%
    median_qi(.width = 0.95) %>%
    select(-c(.width, .point, .interval)) %>%
    mutate(date = recode(date, !!!seq(start_date, now, by = 1))) %>%
    return()
}

summarize_latent <- function(fit, start_date, now, n_lambda_pre, L) {
  fit %>%
    spread_draws(iota[date]) %>%
    median_qi(.width = 0.95) %>%
    select(-c(.width, .point, .interval)) %>%
    mutate(date = recode(date, !!!seq(as.Date(start_date) - n_lambda_pre - L, as.Date(now), by = 1))) %>%
    return()
}

summarize_R <- function(fit, start_date, now, n_lambda_pre, L, max_gen) {
  fit %>%
    spread_draws(R[date]) %>%
    median_qi(.width = 0.95) %>%
    select(-c(.width, .point, .interval)) %>%
    mutate(date = recode(date, !!!seq(as.Date(start_date) - n_lambda_pre - L + max_gen, as.Date(now), by = 1))) %>%
    return()
}

posterior_R <- function(fit, start_date, now, n_lambda_pre, L, max_gen) {
  fit %>%
    spread_draws(R[date]) %>%
    select(-c(.chain, .iteration)) %>%
    mutate(date = recode(date, !!!seq(as.Date(start_date) - n_lambda_pre - L + max_gen, as.Date(now), by = 1))) %>%
    return()
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

draws_R_epiestim <- function(fit,
                             start_date,
                             now,
                             maxT,
                             L,
                             latent_delay_dist,
                             n_samples = 1000,
                             estimation_window = 7,
                             mean_serial_interval = 4.8,
                             std_serial_interval = 2.3,
                             mean_Re_prior = 1,
                             sample_incubation = FALSE) {
  right_bound <- maxT - (estimation_window - 1)
  t_start <- seq(1 + L, right_bound)
  t_end <- t_start + estimation_window - 1

  mean_inc <- weighted.mean(0:L, w = latent_delay_dist)

  regex_args <- quo(`(nowcast_known)|(nowcast_unknown)|(nowcast_all)`[date])

  R_draws <- fit %>%
    spread_draws(!!regex_args, regex = TRUE, ndraws = n_samples)

  if ("nowcast_known" %in% names(R_draws)) {
    R_draws <- R_draws %>% mutate(nowcast_all = nowcast_known)
  }
  if (all(c("nowcast_known", "nowcast_unknown") %in% names(R_draws))) {
    R_draws <- R_draws %>% mutate(nowcast_all = nowcast_known + nowcast_unknown)
  }

  # check if latent delay dist was scaled (to account for reporting proportion != 1)
  rep_prop <- sum(latent_delay_dist) # use further below to scale cases when estimating R
  latent_delay_dist <- latent_delay_dist / rep_prop

  R_draws <- R_draws %>%
    group_by(.draw) %>%
    arrange(date) %>%
    group_modify(~ .sample_R(.x$nowcast_all / rep_prop, mean_serial_interval, std_serial_interval, t_start, t_end, mean_Re_prior)) %>%
    mutate(inc_sample = ifelse(rep(sample_incubation, n()), sample(x = 0:L, size = n(), prob = latent_delay_dist, replace = T), as.integer(ceiling(mean_inc)))) %>%
    transmute(date = recode(t, !!!seq(as.Date(start_date), as.Date(now), by = 1)) - inc_sample, R) %>%
    mutate(date = date - as.integer(floor(estimation_window / 2))) # center R estimates on smoothing window

  return(R_draws)
}

get_imputed_posterior <- function(fitted_model, D, ndraws = 100) {
  df <- fitted_model %>%
    spread_draws(reported_unknown_imputed[date, d], ndraws = ndraws) %>%
    select(.draw, date, d, reported_unknown_imputed) %>%
    mutate(date = date + D) %>%
    arrange(d, date, .draw)

  .draw_unique <- unique(df$.draw)
  date_unique <- unique(df$date)
  d_unique <- unique(df$d)

  post_draws <- array(
    data = df$reported_unknown_imputed,
    dim = c(
      length(.draw_unique),
      length(date_unique),
      length(d_unique)
    ),
    dimnames = list(.draw_unique, date_unique, d_unique)
  )

  return(post_draws)
}
