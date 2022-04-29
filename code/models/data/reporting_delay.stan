/** ---------------------------------------------------------------------------
Data and prior information for reporting delay model
---------------------------------------------------------------------------- */
  // resolution of baseline hazard
  int delay_idx[D+2]; // mapping from delays to baseline hazards

  // baseline hazard gamma_d
  // prior information
  real gamma_mu[n_delays];
  real gamma_sd[n_delays];

  // changepoint model for daily hazard
  int n_beta; // number of changepoints
  row_vector[n_beta] Z[T+D]; // design matrix for changepoint model
  // hyperpriors
  // note: a prior centered at zero assumes that the delay stays constant
  real beta_mu[n_beta];
  real beta_sd[n_beta];

  // reporting day effects / additional covariates
  int n_eta; // number of reporting day / covariate effects
  matrix[n_delays, n_eta] W[T+D]; // the covariate values
  // hyperpriors
  real eta_mu[n_eta];
  real eta_sd[n_eta];

  
  