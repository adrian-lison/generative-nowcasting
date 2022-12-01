/** ---------------------------------------------------------------------------
Parameters for reporting delay model
---------------------------------------------------------------------------- */

  // baseline hazard gamma_d
  vector[n_delays] gamma;
  
  // changepoint model for daily hazard
  array[beta_random] real<lower=0> beta_sd;
  vector<multiplier = (beta_random ? beta_sd[1] : 1)>[n_beta] beta;
 
  // reporting day / additional covariate effects 
  vector[n_eta] eta;