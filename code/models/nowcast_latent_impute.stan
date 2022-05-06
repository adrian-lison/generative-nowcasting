functions  {
  #include functions/helper_functions.stan
}

data {
  int T; // number of days
  int D; // maximum reporting delay
  int n_delays; // number of different baseline hazards
  
  int reported_known[T, D+1]; // reporting triangle
  int reported_unknown[T]; // observations with unknown occurrence date
  
  // include reporting delay data and prior information
  #include data/reporting_delay.stan
  
  // include latent delay data and prior information
  #include data/latent_delay.stan
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
  //real<lower=0> xi_negbinom; // // over-dispersion on the parameter 1 / sqrt(phi) of the negative binomial
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
  //real phi_negbinom = inv_square(xi_negbinom);
  
  // latent event process
  // AR(1) process on log scale
  iota = exp(ar1_process_noncentered_vec(iota_log_start,iota_log_raw,iota_log_sd));
  
  // occurrence process (convolution)
  for(t in 1:(D+T)) {
    lambda[t] = dot_product(latent_delay_dist,iota[(L+t-L):(L+t)]);
  }
  
  // share of events with known occurrence date process
  // AR(1) process on logit scale
  alpha = inv_logit(ar1_process_noncentered_vec(alpha_logit_start,alpha_logit_raw,alpha_logit_sd));
  
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
  //xi_negbinom ~ normal(0., 1.);
  // in Günther et al. 2021, this was either improper,
  // or phi_negbinom ~ inv_gamma(0.01, 0.01);
  
  // random walk prior for log latent events (iota_log)
  iota_log_sd ~ normal(0,0.5) T[0, ]; // truncated normal
  iota_log_start ~ normal(0,12); // starting prior for AR
  iota_log_raw[1:L+D+T] ~ normal(0,1); // non-centered
  
  // random walk prior for share of events with known occurrence date
  alpha_logit_sd ~ normal(0,0.5) T[0, ]; // truncated normal
  alpha_logit_start ~ normal(0,2); // starting prior
  alpha_logit_raw[1:T] ~ normal(0,1); // non-centered

  // Likelihood
  {
  #include model/likelihood_reported.stan
  }
  
}
generated quantities {
  #include generated_quantities/nowcast.stan
}
