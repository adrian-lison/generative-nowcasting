/** ---------------------------------------------------------------------------
Prior sampling statements for reporting delay model
---------------------------------------------------------------------------- */

  // baseline hazard gamma_d
  gamma ~ normal(gamma_mu, gamma_sd);
  
  // changepoint model for daily hazard
  beta ~ normal(beta_mu, beta_sd);

  // reporting day / additional covariate effects
  eta ~ normal(eta_mu, eta_sd);