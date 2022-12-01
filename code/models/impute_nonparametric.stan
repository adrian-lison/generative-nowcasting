functions {
#include functions/helper_functions.stan
#include functions/reporting_triangles.stan  
#include functions/hazards.stan
}

data {
  int T; // number of days
  int D; // maximum reporting delay
  int n_delays; // number of different baseline hazards
  
  int<lower=0> n_lambda_pre; // only used for compatibility with other models
  
  array[T, D+1] int reported_known; // reporting triangle
  array[T] int reported_unknown;
  
  // resolution of baseline hazard
  array[D+2] int delay_idx; // mapping from delays to baseline hazards

  // baseline hazard gamma_d
  // prior information
  array[n_delays] real gamma_prior_mu;
  array[n_delays] real<lower=0> gamma_prior_sd;

  // changepoint model for daily hazard
  int n_beta; // number of changepoints
  array[n_lambda_pre+T] row_vector[n_beta] Z; // design matrix for changepoint model
  // hyperpriors
  // note: a prior centered at zero assumes that the delay stays constant
  array[n_beta] real beta_prior_mu;
  array[n_beta] real<lower=0> beta_prior_sd;
  int<lower=0,upper=1> beta_random;
  array[beta_random] real beta_sd_prior_mu;
  array[beta_random] real<lower=0> beta_sd_prior_sd;

  // reporting day effects / additional covariates
  int n_eta; // number of reporting day / covariate effects
  array[n_lambda_pre+T] matrix[n_delays, n_eta] W; // the covariate values
  // hyperpriors
  array[n_eta] real eta_prior_mu;
  array[n_eta] real<lower=0> eta_prior_sd;
}

transformed data {
  array[T-D, D+1] int reported_known_report = reporting_triangle_by_report(reported_known, D);
}

parameters {
  // baseline hazard gamma_d
  vector[n_delays] gamma;
  // changepoint model for daily hazard
  array[beta_random] real<lower=0> beta_sd;
  vector<multiplier = (beta_random ? beta_sd[1] : 1)>[n_beta] beta;
  // reporting day / additional covariate effects 
  vector[n_eta] eta;
}

transformed parameters {
  
  // delay distribution (order of dimensions reversed to avoid transposition)
  matrix[D+1, T-D] p;
  
  {
    real contrib_occurrence_covariates;
    vector[n_delays] logit_haz_interval;
    vector[D] hazard;
    
    // inference phase
    for(t in (D+1):T) {
      //delay distribution
      contrib_occurrence_covariates = Z[n_lambda_pre+t] * beta;
      logit_haz_interval[1:n_delays] = gamma[1:n_delays] + rep_vector(contrib_occurrence_covariates,n_delays) + W[n_lambda_pre+t,1:n_delays] * eta;
      hazard[1:D] = inv_logit(logit_haz_interval)[delay_idx[1:D]];
  
      p[:, t-D] = compute_prob_from_hazard(hazard[1:D]);
    }
  }

}

model {
  // Priors
  // baseline hazard gamma_d
  gamma ~ normal(gamma_prior_mu, gamma_prior_sd);
  
  // changepoint model for daily hazard
  if (beta_random) {
    beta_sd[1] ~ normal(beta_sd_prior_mu[1], beta_sd_prior_sd[1]) T[0, ]; // truncated normal
    beta ~ normal(beta_prior_mu, beta_sd[1]);
  } else {
    beta ~ normal(beta_prior_mu, beta_prior_sd);
  }

  // reporting day / additional covariate effects
  eta ~ normal(eta_prior_mu, eta_prior_sd);
  
  // Likelihood
  for (t in (D+1):T) {
    reported_known_report[t-D, :] ~ multinomial(p[:, t-D]);
  }
}

generated quantities {
  array[T-D, D+1] int reported_unknown_imputed = rep_array(0, T-D, D+1); // imputed reporting triangle
  for (t in (D+1):T) {
    if (reported_unknown[t] > 0) {
      array[D+1] int backward = multinomial_rng(p[:, t-D], reported_unknown[t]);
      for (d in 1:min(D+1,t-D)) {
        reported_unknown_imputed[t - d + 1 - D, d] += backward[d];
      }
    }
  }
}

