/** ---------------------------------------------------------------------------
Prior sampling statements for reporting delay model
---------------------------------------------------------------------------- */

  // baseline hazard gamma_d
  gamma ~ normal(gamma_prior_mu, gamma_prior_sd);
  
  // changepoint model for daily hazard
  beta ~ normal(beta_prior_mu, beta_prior_sd);

  // reporting day / additional covariate effects
  eta ~ normal(eta_prior_mu, eta_prior_sd);
  
  // random walk prior for share of events with known occurrence date
  alpha_logit_start ~ normal(alpha_logit_start_prior_mu,alpha_logit_start_prior_sd); // starting prior
  alpha_logit_sd ~ normal(alpha_logit_sd_prior_mu,alpha_logit_sd_prior_sd) T[0, ]; // truncated normal
  alpha_logit_noise ~ normal(0, alpha_logit_sd); // Gaussian noise