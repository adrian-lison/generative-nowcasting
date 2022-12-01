/** ---------------------------------------------------------------------------
Data and prior information for reporting delay model
---------------------------------------------------------------------------- */
  // resolution of baseline hazard
  array[D+2] int delay_idx; // mapping from delays to baseline hazards

  // baseline hazard gamma_d
  // prior information
  array[n_delays] real gamma_prior_mu;
  array[n_delays] real<lower=0> gamma_prior_sd;

  // changepoint model for daily hazard
  int n_beta; // number of changepoints
  array[n_lambda_pre+T] row_vector[n_beta] Z; // design matrix for changepoint model
  // hyperpriors
  // note: a prior centered at zero assumes that the delay stays constant
  array[n_beta] real beta_prior_mu;
  array[n_beta] real<lower=0> beta_prior_sd;
  int<lower=0,upper=1> beta_random;
  array[beta_random] real beta_sd_prior_mu;
  array[beta_random] real<lower=0> beta_sd_prior_sd;

  // reporting day effects / additional covariates
  int n_eta; // number of reporting day / covariate effects
  array[n_lambda_pre+T] matrix[n_delays, n_eta] W; // the covariate values
  // hyperpriors
  array[n_eta] real eta_prior_mu;
  array[n_eta] real<lower=0> eta_prior_sd;
  
  // priors for share of events with known occurrence date
  // AR(1) noise sd (truncated normal)
  real alpha_logit_sd_prior_mu; 
  real<lower=0> alpha_logit_sd_prior_sd;
  // AR(1) starting level
  real alpha_logit_start_prior_mu;
  real<lower=0> alpha_logit_start_prior_sd;

  
  