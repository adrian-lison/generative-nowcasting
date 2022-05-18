/** ---------------------------------------------------------------------------
Prior sampling statements for reporting delay model
---------------------------------------------------------------------------- */

  // baseline hazard gamma_d
  gamma ~ normal(gamma_mu, gamma_sd);
  
  // changepoint model for daily hazard
  beta ~ normal(beta_mu, beta_sd);

  // reporting day / additional covariate effects
  eta ~ normal(eta_mu, eta_sd);  // random walk prior for share of events with known occurrence date
  alpha_logit_sd ~ normal(alpha_logit_sd_prior_mu,alpha_logit_sd_prior_sd) T[0, ]; // truncated normal
  alpha_logit_start ~ normal(alpha_logit_start_prior_mu,alpha_logit_start_prior_sd); // starting prior
  alpha_logit_raw[1:T] ~ std_normal(); // non-centered