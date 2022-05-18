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
  
  // include latent delay data
  #include data/latent_delay.stan
  
  // maximum generation time
  int max_gen;

  // delay distribution
  // reversed, such that the probability for a delay of one comes last (zero excluded)
  vector<lower=0>[max_gen] generation_time_dist;
  
  real<lower=0> R_ets_alpha_prior_alpha;
  real<lower=0> R_ets_alpha_prior_beta;
  real<lower=0> R_ets_beta_prior_alpha;
  real<lower=0> R_ets_beta_prior_beta;
  real<lower=0> R_ets_phi_prior_alpha;
  real<lower=0> R_ets_phi_prior_beta;
  real R_sd_prior_mu;
  real<lower=0> R_sd_prior_sd;
  real R_level_start_prior_mu;
  real<lower=0> R_level_start_prior_sd;
  real R_trend_start_prior_mu;
  real<lower=0> R_trend_start_prior_sd;
  
  real<lower=0> iota_initial_prior_mu[max_gen];
  real<lower=0> iota_initial_prior_sd[max_gen];
  
  real xi_negbinom_prior_mu;
  real<lower=0> xi_negbinom_prior_sd;
}

parameters {
  // exponential smoothing / innovations state space process for log(R)
  real<lower=0,upper=1> R_ets_alpha; // smoothing parameter for the level
  real<lower=0,upper=1> R_ets_beta; // smoothing parameter for the trend
  real<lower=0,upper=1> R_ets_phi; // dampening parameter of the trend
  real R_level_start; // starting value of the level
  real R_trend_start; // starting value of the trend
  vector[T+L+D] R_raw; // standardized additive errors
  real<lower=0> R_sd; // standard deviation of additive errors
  
  // random walk parameters for alpha
  real alpha_logit_start;
  vector[T] alpha_logit_raw; // nuisance parameter for non-centered parameterization
  real<lower=0> alpha_logit_sd;
  
  // reporting delay model parameters
  #include parameters/reporting_delay.stan
  
  // initial expected infections
  vector<lower=0>[max_gen] iota_initial;
  
  // realized latent events
  vector<lower=0>[max_gen+L+D+T] I;
  
  // over-dispersion parameter for negative binomial
  real<lower=0> xi_negbinom; // // over-dispersion on the parameter 1 / sqrt(phi) of the negative binomial
}
transformed parameters {
  // expected number of events by occurrence date
  vector<lower=0>[D+T] lambda;
  
  // expected number of latent events by latent date
  vector<lower=0>[L+D+T] iota;
  
  // effective reproduction number
  vector<lower=0>[L+D+T] R;
  
  // delay distribution (order of dimensions reversed to avoid transposition)
  matrix[D+1, T] p;
  matrix[D+1, D] p_pre;
  
  // share of events with known occurrence date, by occurrence date
  vector<lower=0,upper=1>[T] alpha;
  
  // over-dispersion parameter for negative binomial
  real phi_negbinom = inv_square(xi_negbinom);
  
  // Smoothing prior for R
  // ETS/Innovations state space process on log scale, starting value 1 on unit scale
  R = softplus(holt_damped_process_noncentered(R_ets_alpha,R_ets_beta,R_ets_phi,R_level_start,R_trend_start,R_raw,R_sd),4);
  
  // latent event process (convolution) / renewal equation
  for(t in 1:(L+D+T)) {
    iota[t] = R[t] * dot_product(generation_time_dist,I[(max_gen+t-max_gen):(max_gen+t-1)]);
  }
  
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
  xi_negbinom ~ normal(xi_negbinom_prior_mu, xi_negbinom_prior_sd);
  // in Günther et al. 2021, this was either improper,
  // or phi_negbinom ~ inv_gamma(0.01, 0.01);
  
  // ETS/Innovations state space process prior for log R
  R_ets_alpha ~ beta(R_ets_alpha_prior_alpha,R_ets_alpha_prior_beta); 
  R_ets_beta ~ beta(R_ets_beta_prior_alpha,R_ets_beta_prior_beta); 
  R_ets_phi ~ beta(R_ets_phi_prior_alpha,R_ets_phi_prior_beta); // dampening needs a tight prior, roughly between 0.8 and 0.98
  R_sd ~ normal(R_sd_prior_mu,R_sd_prior_sd) T[0, ]; // truncated normal
  R_level_start ~ normal(R_level_start_prior_mu,R_level_start_prior_sd); // starting prior for level
  R_trend_start ~ normal(R_trend_start_prior_mu,R_trend_start_prior_sd); // starting prior for trend
  R_raw[1:L+D+T] ~ std_normal(); // non-centered
  
  // latent event realizations
  iota_initial ~ normal(iota_initial_prior_mu,iota_initial_prior_sd); // half normal due to constraint
  I[1:max_gen] ~ normal(iota_initial,sqrt(iota_initial)); // half normal due to constraint, approximates Poisson
  I[(max_gen+1):(max_gen+L+D+T)] ~ normal(iota,sqrt(iota)); // half normal due to constraint, approximates Poisson

  // Likelihood
  {
  #include model/likelihood_reported.stan
  }
  
}
generated quantities {
  #include generated_quantities/nowcast.stan
}
