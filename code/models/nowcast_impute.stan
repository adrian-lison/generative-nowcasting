functions  {
#include functions/helper_functions.stan
#include functions/time_series.stan
#include functions/hazards.stan
#include functions/reporting_triangles.stan
#include functions/distribution_approximations.stan
}

data {
  int T; // number of days
  int D; // maximum reporting delay
  int n_delays; // number of different baseline hazards
  
  int<lower=0> n_lambda_pre; // number of days used for modeling of latent process before first data
  
  array[T, D+1] int reported_known; // reporting triangle
  array[T] int reported_unknown; // observations with unknown occurrence date
  
  // include reporting delay data and prior information
#include data/reporting_delay.stan
  
  // include priors/settings for exponential smoothing
#include data/data_ets.stan
  
  // standard deviation of noise/increments
  real lambda_log_sd_prior_mu;
  real<lower=0> lambda_log_sd_prior_sd;
  // intercept
  real lambda_log_level_start_prior_mu;
  real<lower=0> lambda_log_level_start_prior_sd;
  // trend
  real lambda_log_trend_start_prior_mu;
  real<lower=0> lambda_log_trend_start_prior_sd;
  // second order trend / curvature
  array[ets_diff ? 1 : 0] real lambda_log_2nd_trend_start_prior_mu;
  array[ets_diff ? 1 : 0] real<lower=0> lambda_log_2nd_trend_start_prior_sd;
  
  int<lower=0,upper=1> overdispersion; // whether to model overdispersion via a negative binomial
  real xi_negbinom_prior_mu;
  real<lower=0> xi_negbinom_prior_sd;
}

parameters {
  // exponential smoothing / innovations state space process for log(lambda)
  array[ets_alpha_fixed < 0 ? 1 : 0] real<lower=0,upper=1> ets_alpha; // smoothing parameter for the level
  array[ets_beta_fixed < 0 ? 1 : 0] real<lower=0,upper=1> ets_beta; // smoothing parameter for the trend
  array[ets_phi_fixed < 0 ? 1 : 0] real<lower=0,upper=1> ets_phi; // dampening parameter of the trend
  vector[2+ets_diff] lambda_log_start_values; // starting value of the level, trend [, and curvature if diff]
  real<lower=0> lambda_log_sd; // standard deviation of process noise
  vector<multiplier=(ets_noncentered ? lambda_log_sd : 1)>[n_lambda_pre+T-1-ets_diff] lambda_log_noise; // process noise
  
  // random walk parameters for alpha
  real alpha_logit_start;
  real<lower=0> alpha_logit_sd;
  vector<multiplier=alpha_logit_sd>[T-1] alpha_logit_noise; // process noise
  
  // reporting delay model parameters
#include parameters/reporting_delay.stan
  
  // over-dispersion on the parameter 1 / sqrt(phi) of the negative binomial
  array[overdispersion ? 1 : 0] real<lower=0> xi_negbinom;
}
transformed parameters {
  // expected number of events by occurrence date
  vector[n_lambda_pre+T] lambda_log; 
  
  // delay distribution (order of dimensions reversed to avoid transposition)
  matrix[D+1, T] p_log;
  matrix[D+1, n_lambda_pre] p_log_pre;
  
  // share of events with known occurrence date, by occurrence date
  vector[T] alpha_log;
  vector[T] alpha1m_log;
  
  // over-dispersion parameter for negative binomial
  real phi_negbinom;
  if (overdispersion) {
    phi_negbinom = inv_square(xi_negbinom[1]);
  }
  
  // occurrence process
  lambda_log = holt_damped_process(
    lambda_log_start_values,
    ets_alpha_fixed < 0 ? ets_alpha[1] : ets_alpha_fixed,
    ets_beta_fixed < 0 ? ets_beta[1] : ets_beta_fixed,
    ets_phi_fixed < 0 ? ets_phi[1] : ets_phi_fixed,
    lambda_log_noise,
    ets_diff);
    
  // share of events with known occurrence date process
  // AR(1) process on logit scale
  {
    vector[T] alpha_logit = random_walk([alpha_logit_start]', alpha_logit_noise, 0);
    alpha_log = log_inv_logit(alpha_logit);
    alpha1m_log = log1m_inv_logit(alpha_logit);
  }

  
  // reporting delay model: hazards and probabilities
  {
#include transformed_parameters/reporting_delay.stan
  }
}

model {
  // Priors
  
  // priors for reporting delay model
#include model/priors_reporting_delay.stan
  
  // prior for overdispersion
  if (overdispersion) {
    xi_negbinom[1] ~ normal(xi_negbinom_prior_mu, xi_negbinom_prior_sd) T[0, ]; // truncated normal
  }
  // in Günther et al. 2021, this was either improper,
  // or phi_negbinom ~ inv_gamma(0.01, 0.01);
  
  // ETS/Innovations state space process prior for log(lambda)
  if(ets_alpha_fixed < 0) {
    ets_alpha[1] ~ beta(ets_alpha_prior_alpha[1], ets_alpha_prior_beta[1]); 
  }
  if(ets_beta_fixed < 0) {
    ets_beta[1] ~ beta(ets_beta_prior_alpha[1], ets_beta_prior_beta[1]); 
  }
  if(ets_phi_fixed < 0) {
    ets_phi[1] ~ beta(ets_phi_prior_alpha[1], ets_phi_prior_beta[1]); // dampening needs a tight prior, roughly between 0.8 and 0.98
  }
  lambda_log_start_values[1] ~ normal(lambda_log_level_start_prior_mu, lambda_log_level_start_prior_sd); // starting prior for AR
  lambda_log_start_values[2] ~ normal(lambda_log_trend_start_prior_mu, lambda_log_trend_start_prior_sd);
  if (ets_diff == 1) {
    lambda_log_start_values[3] ~ normal(lambda_log_2nd_trend_start_prior_mu, lambda_log_2nd_trend_start_prior_sd);
  }
  lambda_log_sd ~ normal(lambda_log_sd_prior_mu, lambda_log_sd_prior_sd) T[0, ]; // truncated normal
  lambda_log_noise ~ normal(0, lambda_log_sd); // Gaussian noise

  // Likelihood
  {
#include model/likelihood_reported.stan
  }
  
}
generated quantities {
#include generated_quantities/nowcast.stan
}
