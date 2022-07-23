/** ---------------------------------------------------------------------------
Data and prior information for reporting delay model
---------------------------------------------------------------------------- */
  // resolution of baseline hazard
  int delay_idx[D+2]; // mapping from delays to baseline hazards

  // baseline hazard gamma_d
  // prior information
  real gamma_prior_mu[n_delays];
  real<lower=0> gamma_prior_sd[n_delays];

  // changepoint model for daily hazard
  int n_beta; // number of changepoints
  row_vector[n_beta] Z[n_lambda_pre+T]; // design matrix for changepoint model
  // hyperpriors
  // note: a prior centered at zero assumes that the delay stays constant
  real beta_prior_mu[n_beta];
  real<lower=0> beta_prior_sd[n_beta];

  // reporting day effects / additional covariates
  int n_eta; // number of reporting day / covariate effects
  matrix[n_delays, n_eta] W[n_lambda_pre+T]; // the covariate values
  // hyperpriors
  real eta_prior_mu[n_eta];
  real<lower=0> eta_prior_sd[n_eta];
  
  // priors for share of events with known occurrence date
  // AR(1) noise sd (truncated normal)
  real alpha_logit_sd_prior_mu; 
  real<lower=0> alpha_logit_sd_prior_sd;
  // AR(1) starting level
  real alpha_logit_start_prior_mu;
  real<lower=0> alpha_logit_start_prior_sd;

  
  