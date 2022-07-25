/** ---------------------------------------------------------------------------
Priors / settings for exponential smoothing / latent state space model
---------------------------------------------------------------------------- */
  int<lower=0> ets_diff; // order of differencing
  int<lower=0, upper=1> ets_noncentered; // use non-centered parameterization?
  
  real ets_alpha_fixed; // fixed value used if non-negative
  array[ets_alpha_fixed < 0 ? 1 : 0] real<lower=0> ets_alpha_prior_alpha;
  array[ets_alpha_fixed < 0 ? 1 : 0] real<lower=0> ets_alpha_prior_beta;
  
  real ets_beta_fixed; // fixed value used if non-negative
  array[ets_beta_fixed < 0 ? 1 : 0] real<lower=0> ets_beta_prior_alpha;
  array[ets_beta_fixed < 0 ? 1 : 0] real<lower=0> ets_beta_prior_beta;
  
  real ets_phi_fixed; // fixed value used if non-negative
  array[ets_phi_fixed < 0 ? 1 : 0] real<lower=0> ets_phi_prior_alpha;
  array[ets_phi_fixed < 0 ? 1 : 0] real<lower=0> ets_phi_prior_beta;