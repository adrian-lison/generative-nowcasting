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
  
  int n_imputations; // number of different imputations
  array[n_imputations, T, D+1] int reported_unknown_imputed; // imputed reporting triangle
  
  // include reporting delay data and prior information
#include data/reporting_delay.stan

  // include latent delay data
#include data/latent_delay.stan
  
  // maximum generation time
  int max_gen;

  // delay distribution
  // probability for a delay of one comes first (zero excluded)
  vector<lower=0>[max_gen] generation_time_dist;
  
  // include priors/settings for exponential smoothing
#include data/data_ets.stan
  int<lower=0> ts_model; 
  int<lower=0> sma_window;
  
  real R_sd_prior_mu;
  real<lower=0> R_sd_prior_sd;
  real R_level_start_prior_mu;
  real<lower=0> R_level_start_prior_sd;
  real R_trend_start_prior_mu;
  real<lower=0> R_trend_start_prior_sd;
  
  real iota_log_ar_start_prior_mu;
  real<lower=0> iota_log_ar_start_prior_sd;
  real iota_log_ar_sd_prior_mu;
  real<lower=0> iota_log_ar_sd_prior_sd;
  
  int<lower=0,upper=1> overdispersion; // whether to model overdispersion via a negative binomial
  real xi_negbinom_prior_mu;
  real<lower=0> xi_negbinom_prior_sd;
}

transformed data {
  vector[max_gen] generation_time_dist_reversed;
  vector[L+1] latent_delay_dist_reversed;
  array[n_imputations, T, D+1] int reported_imputed;
  
  generation_time_dist_reversed = reverse(generation_time_dist);
  latent_delay_dist_reversed = reverse(latent_delay_dist);
  
  for (i in 1:n_imputations) {
    for (t in 1:T) {
      for (d in 1:(D+1)) {
        reported_imputed[i,t,d] = reported_known[t,d] + reported_unknown_imputed[i,t,d];
      }
    }
  }
}

parameters {
  /// log(R) time series prior
  real R_level_start; // starting value of the level
  real R_trend_start; // starting value of the trend
  real<lower=0> R_sd; // standard deviation of additive errors
  vector<multiplier=(ts_model > 0 ? (ets_noncentered ? R_sd : 1) : 1)>[L+n_lambda_pre+T-max_gen-1] R_noise; // additive errors
  
  // exponential smoothing / innovations state space process for log(R)
  array[ts_model == 1 ? (ets_alpha_fixed < 0 ? 1 : 0) : 0] real<lower=0,upper=1> ets_alpha; // smoothing parameter for the level
  array[ts_model == 1 ? (ets_beta_fixed < 0 ? 1 : 0) : 0] real<lower=0,upper=1> ets_beta; // smoothing parameter for the trend
  array[ts_model == 1 ? (ets_phi_fixed < 0 ? 1 : 0) : 0] real<lower=0,upper=1> ets_phi; // dampening parameter of the trend
  
  // reporting delay model parameters
#include parameters/reporting_delay.stan

  // realized latent events
  real iota_log_ar_start;
  real<lower=0> iota_log_ar_sd;
  vector<multiplier=iota_log_ar_sd>[max_gen-1] iota_log_ar_noise;
  vector<lower=0>[L+n_lambda_pre+T] I;
  
  // over-dispersion on the parameter 1 / sqrt(phi) of the negative binomial
  array[overdispersion ? 1 : 0] real<lower=0> xi_negbinom;
}
transformed parameters {
  // expected number of events by occurrence date
  vector[n_lambda_pre+T] lambda_log; 
  
  // expected number of latent events by latent date
  vector[L+n_lambda_pre+T] iota;
  
  // effective reproduction number
  vector[L+n_lambda_pre+T-max_gen] R;
  
  // delay distribution (order of dimensions reversed to avoid transposition)
  matrix[D+1, T] p_log;
  matrix[D+1, n_lambda_pre] p_log_pre;
  
  // over-dispersion parameter for negative binomial
  real phi_negbinom;
  if (overdispersion) {
    phi_negbinom = inv_square(xi_negbinom[1]);
  }
  
  // Smoothing prior for R
  // ETS/Innovations state space process on log scale, starting value 1 on unit scale
  profile ("transformed_R") {
  if (ts_model == 0) {
    R = softplus(append_row(R_level_start, R_noise), 4);
  } else if (ts_model == 1) {
    R = softplus(holt_damped_process(
      [R_level_start, R_trend_start]',
      ets_alpha_fixed < 0 ? ets_alpha[1] : ets_alpha_fixed,
      ets_beta_fixed < 0 ? ets_beta[1] : ets_beta_fixed,
      ets_phi_fixed < 0 ? ets_phi[1] : ets_phi_fixed,
      R_noise, 0), 4);
  } else if (ts_model == 2) {
    R = softplus(simple_ma([R_level_start, R_trend_start]', R_noise, sma_window, 0), 4);
  }
  }
  
  // latent event process (convolution) / renewal equation
  profile ("transformed_infections") {
  // for t in 1:max_gen
  // vector[L+n_lambda_pre+T] I_log = log(I+0.01);
  iota[1:max_gen] = exp(random_walk([iota_log_ar_start]', iota_log_ar_noise, 0));
  for(t in (max_gen+1):(L+n_lambda_pre+T)) {
    iota[t] = R[t-max_gen] * dot_product(generation_time_dist_reversed,I[(t-max_gen):(t-1)]);
  }
  }
  
  // occurrence process (convolution)
  profile ("transformed_lambda") {
  for(t in 1:(n_lambda_pre+T)) {
    lambda_log[t] = log(dot_product(latent_delay_dist_reversed,iota[(L+t-L):(L+t)]));
  }
  }
  
  profile ("transformed_reporting") {
  // reporting delay model: hazards and probabilities
  {
#include transformed_parameters/reporting_delay.stan
  }
  }
}

model {
  // Priors
  
    // baseline hazard gamma_d
  gamma ~ normal(gamma_prior_mu, gamma_prior_sd);
  // changepoint model for daily hazard
  beta ~ normal(beta_prior_mu, beta_prior_sd);
  // reporting day / additional covariate effects
  eta ~ normal(eta_prior_mu, eta_prior_sd);
  
  // prior for overdispersion
  if (overdispersion) {
    xi_negbinom[1] ~ normal(xi_negbinom_prior_mu, xi_negbinom_prior_sd) T[0, ]; // truncated normal
  }
  
  // ETS/Innovations state space process prior for log R
  if (ts_model == 1) {
    if(ets_alpha_fixed < 0) {
      ets_alpha[1] ~ beta(ets_alpha_prior_alpha[1], ets_alpha_prior_beta[1]); 
    }
    if(ets_beta_fixed < 0) {
      ets_beta[1] ~ beta(ets_beta_prior_alpha[1], ets_beta_prior_beta[1]); 
    }
    if(ets_phi_fixed < 0) {
      ets_phi[1] ~ beta(ets_phi_prior_alpha[1], ets_phi_prior_beta[1]); // dampening needs a tight prior, roughly between 0.8 and 0.98
    }
  }
  if (ts_model == 0) {
    R_level_start ~ normal(R_level_start_prior_mu, R_level_start_prior_sd); // starting prior for level
    R_noise ~ normal(R_level_start_prior_mu, R_level_start_prior_sd); // independent sampling
  } else {
    R_level_start ~ normal(R_level_start_prior_mu, R_level_start_prior_sd); // starting prior for level
    R_trend_start ~ normal(R_trend_start_prior_mu, R_trend_start_prior_sd); // starting prior for trend
    R_sd ~ normal(R_sd_prior_mu, R_sd_prior_sd) T[0, ]; // truncated normal
    R_noise ~ normal(0, R_sd); // Gaussian noise
  }
  
  // latent event realizations
  iota_log_ar_start ~ normal(iota_log_ar_start_prior_mu, iota_log_ar_start_prior_sd);
  iota_log_ar_sd ~ normal(iota_log_ar_sd_prior_mu, iota_log_ar_sd_prior_sd) T[0, ]; // truncated normal
  iota_log_ar_noise ~ normal(0, iota_log_ar_sd); // Gaussian noise
  I[1:(L+n_lambda_pre+T)] ~ normal(iota, sqrt(iota)); // half normal due to constraint, approximates Poisson

  // Likelihood
  {
#include model/likelihood_reported_multiple_impute.stan
  }
  
}
generated quantities {
#include generated_quantities/nowcast_multiple_impute.stan
}
