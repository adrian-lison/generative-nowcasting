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
  real xi_negbinom_prior_mu;
  real<lower=0> xi_negbinom_prior_sd;
}

parameters {
  // random walk parameters for lambda
  real lambda_log_start;
  vector[D+T] lambda_log_raw; // nuisance parameter for non-centered parameterization
  real<lower=0> lambda_log_sd;
  
  // random walk parameters for alpha
  real alpha_logit_start;
  vector[T] alpha_logit_raw; // nuisance parameter for non-centered parameterization
  real<lower=0> alpha_logit_sd;
  
  // reporting delay model parameters
  #include parameters/reporting_delay.stan
  
  // over-dispersion on the parameter 1 / sqrt(phi) of the negative binomial
  real<lower=0> xi_negbinom;
}
transformed parameters {
  // expected number of events by occurrence date
  vector[D+T] lambda; 
  
  // delay distribution (order of dimensions reversed to avoid transposition)
  matrix[D+1, T] p;
  matrix[D+1, D] p_pre;
  
  // share of events with known occurrence date, by occurrence date
  vector<lower=0,upper=1>[T] alpha;
  
  // over-dispersion parameter for negative binomial
  real phi_negbinom = inv_square(xi_negbinom);
  
  // occurrence process
  // AR(1) process on log scale
  lambda = exp(ar1_process_noncentered_vec(lambda_log_start,lambda_log_raw,lambda_log_sd));
  
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
  xi_negbinom ~ normal(xi_negbinom_prior_mu, xi_negbinom_prior_sd);
  // in Günther et al. 2021, this was either improper,
  // or phi_negbinom ~ inv_gamma(0.01, 0.01);
  
  // random walk prior for log(lambda)
  lambda_log_sd ~ normal(0,0.5) T[0, ]; // truncated normal
  lambda_log_start ~ normal(0,12); // starting prior for AR
  lambda_log_raw[1:(D+T)] ~ normal(0,1); // non-centered


  // Likelihood
  {
  #include model/likelihood_reported.stan
  }
  
}
generated quantities {
  #include generated_quantities/nowcast.stan
}
