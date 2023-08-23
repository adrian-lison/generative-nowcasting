functions {
  array[] int negative_binomial_2_ub_n_rng(int n, real mu, real phi, int bound) {
    array[n] int y = rep_array(0, n);
    vector[bound+1] negbincdf;
    for (i in 0:bound) {
      negbincdf[i+1] = neg_binomial_2_cdf(i | mu, phi);
    }
    real p_ub = negbincdf[bound+1];
    for (i in 1:n) {
      real u = uniform_rng(0, p_ub);
      // inverse cdf
      for (j in 0:bound) {
        if (u > negbincdf[j + 1]) {
          y[i] += 1;
        } else {
          break;
        }
      }
    }
    return y;
  }
}

data {
  int T; // number of days
  int D; // maximum reporting delay
  int n_delays; // number of different baseline hazards
  
  array[T, D+1] int reported_known; // reporting triangle
  array[T] int reported_unknown;
  
  real intercept_prior_mu;
  real<lower=0> intercept_prior_sd;
  
  // changepoint model for daily hazard
  int n_beta; // number of changepoints
  matrix[T, n_beta] Z; // design matrix for changepoint model
  // hyperpriors
  // note: a prior centered at zero assumes that the delay stays constant
  vector[n_beta] beta_prior_mu;
  vector<lower=0>[n_beta] beta_prior_sd;

  // reporting day effects / additional covariates
  int n_eta; // number of reporting day / covariate effects
  array[T] matrix[n_delays, n_eta] W; // the covariate values
  // hyperpriors
  vector[n_eta] eta_prior_mu;
  vector<lower=0>[n_eta] eta_prior_sd;
  
  real xi_negbinom_prior_mu;
  real<lower=0> xi_negbinom_prior_sd;
}

transformed data {
  matrix[T, n_eta] W_mat = to_matrix(W[,1,]);
}

parameters {
  real intercept;
  
  // changepoint model for daily hazard
  vector[n_beta] beta;
 
  // reporting day / additional covariate effects 
  vector[n_eta] eta;
  
  // over-dispersion parameter for negative binomial
  real<lower=0> xi_negbinom;
}

transformed parameters {
  vector[T-D] mu = exp(intercept + Z[(D+1):T,] * beta + W_mat[(D+1):T,] * eta);
  real phi_negbinom = inv_square(xi_negbinom);
}

model {
  // Priors
  intercept ~ normal(intercept_prior_mu, intercept_prior_sd);
  
  // changepoint model for daily hazard
  beta ~ normal(beta_prior_mu, beta_prior_sd);

  // reporting day / additional covariate effects
  eta ~ normal(eta_prior_mu, eta_prior_sd);
  
  // overdispersion
  xi_negbinom ~ normal(xi_negbinom_prior_mu, xi_negbinom_prior_sd) T[0, ]; // truncated normal
  
  // Likelihood
  real d_loglik;
  for (t in 1:T) {
    for (d in max(0,D-t+1):min(D, T-t)) {
      d_loglik = neg_binomial_2_lpmf(d | mu[t - D + d], phi_negbinom);
      d_loglik += - neg_binomial_2_lcdf(D | mu[t - D + d], phi_negbinom); // adjust for right-truncation
      target += reported_known[t, d + 1] * d_loglik; // multiply by number of observations with that delay
    }
  }
}

generated quantities {
  array[T-D, D+1] int reported_unknown_imputed = rep_array(0, T-D, D+1); // imputed reporting triangle
  for (t in (D+1):T) {
    int n = reported_unknown[t];
    if (n>0) {
     array[n] int delays = negative_binomial_2_ub_n_rng(n, mu[t - D], phi_negbinom, D);
    for (i in 1:n) {
      int ref_time = t - delays[i];
      if (ref_time > D) {
        reported_unknown_imputed[ref_time - D, delays[i] + 1] += 1;
      }
    }
    } 
    }
}

