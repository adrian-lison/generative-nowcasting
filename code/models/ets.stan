functions {
  #include functions/helper_functions.stan
  #include functions/time_series.stan
  #include functions/hazards.stan
  #include functions/reporting_triangles.stan
  #include functions/distribution_approximations.stan
}

data {
  int T; // number of observations
  real<lower=0> obs[T]; // observations
  
  real R_ets_alpha_fixed; // fixed value used if non-negative
  real<lower=0> R_ets_alpha_prior_alpha[R_ets_alpha_fixed < 0 ? 1 : 0];
  real<lower=0> R_ets_alpha_prior_beta[R_ets_alpha_fixed < 0 ? 1 : 0];
  real R_ets_beta_fixed; // fixed value used if non-negative
  real<lower=0> R_ets_beta_prior_alpha[R_ets_beta_fixed < 0 ? 1 : 0];
  real<lower=0> R_ets_beta_prior_beta[R_ets_beta_fixed < 0 ? 1 : 0];
  real R_ets_phi_fixed; // fixed value used if non-negative
  real<lower=0> R_ets_phi_prior_alpha[R_ets_phi_fixed < 0 ? 1 : 0];
  real<lower=0> R_ets_phi_prior_beta[R_ets_phi_fixed < 0 ? 1 : 0];
  real R_sd_prior_mu;
  real<lower=0> R_sd_prior_sd;
  real R_level_start_prior_mu;
  real<lower=0> R_level_start_prior_sd;
  real R_trend_start_prior_mu;
  real<lower=0> R_trend_start_prior_sd;
}

parameters {
  // exponential smoothing / innovations state space process for R
  real<lower=0,upper=1> R_ets_alpha[R_ets_alpha_fixed < 0 ? 1 : 0]; // smoothing parameter for the level
  real<lower=0,upper=1> R_ets_beta[R_ets_beta_fixed < 0 ? 1 : 0]; // smoothing parameter for the trend
  real<lower=0,upper=1> R_ets_phi[R_ets_phi_fixed < 0 ? 1 : 0]; // dampening parameter of the trend
  real R_level_start; // starting value of the level
  real R_trend_start; // starting value of the trend
  vector[T] R_raw; // standardized additive errors
  real<lower=0> R_sd; // standard deviation of additive errors
}
transformed parameters {
    // effective reproduction number
  vector<lower=0>[T] R;

  // ETS/Innovations state space process
  R = softplus(holt_damped_process_noncentered(
    R_ets_alpha_fixed < 0 ? R_ets_alpha[1] : R_ets_alpha_fixed,
    R_ets_beta_fixed < 0 ? R_ets_beta[1] : R_ets_beta_fixed,
    R_ets_phi_fixed < 0 ? R_ets_phi[1] : R_ets_phi_fixed,
    R_level_start,R_trend_start,R_raw,R_sd),4);
}

model {
  // Priors
  // ETS/Innovations state space process prior
  if(R_ets_alpha_fixed < 0) {
    R_ets_alpha[1] ~ beta(R_ets_alpha_prior_alpha[1],R_ets_alpha_prior_beta[1]); 
  }
  if(R_ets_beta_fixed < 0) {
    R_ets_beta[1] ~ beta(R_ets_beta_prior_alpha[1],R_ets_beta_prior_beta[1]); 
  }
  if(R_ets_phi_fixed < 0) {
    R_ets_phi[1] ~ beta(R_ets_phi_prior_alpha[1],R_ets_phi_prior_beta[1]); // dampening needs a tight prior, roughly between 0.8 and 0.98
  }
  R_sd ~ normal(R_sd_prior_mu,R_sd_prior_sd) T[0, ]; // truncated normal
  R_level_start ~ normal(R_level_start_prior_mu,R_level_start_prior_sd); // starting prior for level
  R_trend_start ~ normal(R_trend_start_prior_mu,R_trend_start_prior_sd); // starting prior for trend
  R_raw ~ std_normal(); // non-centered

  // Likelihood
  obs ~ normal(R, 0.1);
}
