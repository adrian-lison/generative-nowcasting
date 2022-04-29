/** ---------------------------------------------------------------------------
Parameters for reporting delay model
---------------------------------------------------------------------------- */

  // baseline hazard gamma_d
  vector[n_delays] gamma;
  
  // changepoint model for daily hazard
  vector[n_beta] beta;
 
  // reporting day / additional covariate effects 
  vector[n_eta] eta;