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
  
  int reported_known[T, D+1]; // reporting triangle
  int reported_unknown[T]; // observations with unknown occurrence date
  
  // include reporting delay data and prior information
  #include data/reporting_delay.stan
  
  // include latent delay data
  #include data/latent_delay.stan
  
  // priors for iota random walk
  real iota_log_sd_prior_mu;
  real<lower=0> iota_log_sd_prior_sd;
  real iota_log_start_prior_mu;
  real<lower=0> iota_log_start_prior_sd;
  
  real xi_negbinom_prior_mu;
  real<lower=0> xi_negbinom_prior_sd;
}

transformed data {
  vector[L+1] latent_delay_dist_reversed;
  
  latent_delay_dist_reversed = reverse(latent_delay_dist);
}

parameters {
  // random walk parameters for lambda
  real iota_log_start;
  vector[T+L+D] iota_log_raw; // nuisance parameter for non-centered parameterization
  real<lower=0> iota_log_sd;
  
  // random walk parameters for alpha
  real alpha_logit_start;
  vector[T] alpha_logit_raw; // nuisance parameter for non-centered parameterization
  real<lower=0> alpha_logit_sd;
  
  // reporting delay model parameters
  #include parameters/reporting_delay.stan
  
  // over-dispersion parameter for negative binomial
  real xi_negbinom; // // over-dispersion on the parameter 1 / sqrt(phi) of the negative binomial
}
transformed parameters {
  // expected number of events by occurrence date
  vector[D+T] lambda;
  
  // expected number of latent events by latent date
  vector[L+D+T] iota;
  
  // delay distribution (order of dimensions reversed to avoid transposition)
  matrix[D+1, T] p;
  matrix[D+1, D] p_pre;
  
  // share of events with known occurrence date, by occurrence date
  vector<lower=0,upper=1>[T] alpha;
  
  // over-dispersion parameter for negative binomial
  real phi_negbinom = inv_square(xi_negbinom);
  
  // latent event process
  // AR(1) process on log scale
  iota = exp(random_walk(iota_log_start,iota_log_raw,iota_log_sd));
  
  // occurrence process (convolution)
  for(t in 1:(D+T)) {
    lambda[t] = dot_product(latent_delay_dist_reversed,iota[(L+t-L):(L+t)]);
  }
  
  // share of events with known occurrence date process
  // AR(1) process on logit scale
  alpha = inv_logit(random_walk(alpha_logit_start,alpha_logit_raw,alpha_logit_sd));
  
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
  xi_negbinom ~ normal(xi_negbinom_prior_mu, xi_negbinom_prior_sd);
  // in Günther et al. 2021, this was either improper,
  // or phi_negbinom ~ inv_gamma(0.01, 0.01);
  
  // random walk prior for log latent events (iota_log)
  iota_log_sd ~ normal(iota_log_sd_prior_mu,iota_log_sd_prior_sd) T[0, ]; // truncated normal
  iota_log_start ~ normal(iota_log_start_prior_mu,iota_log_start_prior_sd); // starting prior for AR
  iota_log_raw[1:L+D+T] ~ std_normal(); // non-centered

  // Likelihood
  {
  #include model/likelihood_reported.stan
  }
  
}
generated quantities {
  #include generated_quantities/nowcast.stan
}
