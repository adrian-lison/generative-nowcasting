data {
  // Number of days
  int T;
  // Number of categories
  int J;
  // D
  int maxDelay;
  // Number cps delexp(5ay dist
  int n_cp;
  // Number additional Covariate effects delay distribution
  int n_wextra;
  // observation triangle
  array[J, T, maxDelay + 1] int rT;
  // observations with unknown onset date
  array[J, T] int aT;
  // Design matrix delay dist
  // Note: the way the changepoints are parameterized, a prior centered at zero assumes that
  // the linear TREND in delay continues, not that the delay stays constant
  array[T] matrix[maxDelay + 1, T] W_cp;
  array[T] int W_cp_idx; // index of changepoints
  array[T + maxDelay] matrix[maxDelay + 1, n_wextra] W_extra;
  // Lambda Random Walk Params
  // Prior Delay gammas (log-OR Haz)
  array[maxDelay] real mu_gamma;
  array[maxDelay] real sd_gamma;
  // Prior Effects Delay dist (Changepoints and Wextra)
  array[n_cp] real eta_cp_mu;
  array[n_cp] real eta_cp_sd;
  array[n_wextra] real eta_extra_mu;
  array[n_wextra] real eta_extra_sd;
  // Predictions up to predLag days before now
  int predLag;
}
parameters {
  array[J, T] real logLambda_raw; // nuisance parameter for non-centered parameterization
  array[J, maxDelay] real logLambda_pre_raw; // nuisance parameter for non-centered parameterization
  array[J, T] real logitAlpha_raw; // nuisance parameter for non-centered parameterization
  array[maxDelay] real gamma;
  vector[n_cp] eta_cp;
  vector[n_wextra] eta_extra;
  array[J] real<lower=0> sd_logLambda;
  array[J] real<lower=0> sd_logitAlpha;
  //real<lower=0> xi_negbinom; // // over-dispersion on the parameter 1 / sqrt(phi) of the negative binomial
  array[J - 1] real theta; // logit difference in hazard for categories (currently assuming this stays constant over time)
}
transformed parameters {
  array[J, T] real logLambda;
  array[J, maxDelay] real logLambda_pre;
  array[J, T] real logitAlpha;
  array[J, T, maxDelay + 1] real p;
  array[J, maxDelay, maxDelay + 1] real p_pre;
  array[J, T, maxDelay] real haz;
  array[J, maxDelay, maxDelay] real haz_pre;
  array[J, T] real<lower=0, upper=1> alpha;
  
  //real phi_negbinom = inv_square(xi_negbinom); // over-dispersion parameter for negative binomial
  
  for (j in 1 : J) {
    //for each case category j
    
    // AR processes
    logLambda_pre[j, 1] = logLambda_pre_raw[j, 1];
    for (t_pre in 2 : maxDelay) {
      logLambda_pre[j, t_pre] = logLambda_pre[j, t_pre - 1]
                                + sd_logLambda[j]
                                  * logLambda_pre_raw[j, t_pre];
    }
    logLambda[j, 1] = logLambda_pre[j, maxDelay]
                      + sd_logLambda[j] * logLambda_raw[j, 1];
    logitAlpha[j, 1] = logitAlpha_raw[j, 1];
    for (t in 2 : T) {
      logLambda[j, t] = logLambda[j, t - 1]
                        + sd_logLambda[j] * logLambda_raw[j, t];
      logitAlpha[j, t] = logitAlpha[j, t - 1]
                         + sd_logitAlpha[j] * logitAlpha_raw[j, t];
    }
    
    // delay model
    // modeling phase
    for (t_pre in 1 : maxDelay) {
      //delay distribution
      haz_pre[j, t_pre, 1] = inv_logit(gamma[1]
                                       + W_extra[t_pre, 1] * eta_extra
                                       + (j == 1 ? 0 : theta[j - 1]));
      p_pre[j, t_pre, 1] = haz_pre[j, t_pre, 1];
      for (d in 1 : (maxDelay - 1)) {
        haz_pre[j, t_pre, d + 1] = inv_logit(gamma[d + 1]
                                             + W_extra[t_pre, d + 1]
                                               * eta_extra
                                             + (j == 1 ? 0 : theta[j - 1]));
        p_pre[j, t_pre, d + 1] = (1 - sum(p_pre[j, t_pre, 1 : d]))
                                 * haz_pre[j, t_pre, d + 1];
      }
      p_pre[j, t_pre, maxDelay + 1] = 1 - sum(p_pre[j, t_pre, 1 : maxDelay]);
    }
    // inference phase
    for (t in 1 : T) {
      //delay distribution
      haz[j, t, 1] = inv_logit(gamma[1] + W_cp[t, 1] * eta_cp[W_cp_idx]
                               + W_extra[t + maxDelay, 1] * eta_extra
                               + (j == 1 ? 0 : theta[j - 1]));
      p[j, t, 1] = haz[j, t, 1];
      for (d in 1 : (maxDelay - 1)) {
        haz[j, t, d + 1] = inv_logit(gamma[d + 1]
                                     + W_cp[t, d + 1] * eta_cp[W_cp_idx]
                                     + W_extra[t + maxDelay, d + 1]
                                       * eta_extra
                                     + (j == 1 ? 0 : theta[j - 1]));
        p[j, t, d + 1] = (1 - sum(p[j, t, 1 : d])) * haz[j, t, d + 1];
      }
      p[j, t, maxDelay + 1] = 1 - sum(p[j, t, 1 : maxDelay]);
    }
  }
  
  alpha = inv_logit(logitAlpha);
}
model {
  // priors
  // logLambda
  sd_logLambda[ : ] ~ normal(0, 0.5); // half normal due to constraint: made much smaller than original N(0,5), since it is geometric random walk - a log increase of 1 already equals almost a trippling per day
  
  logLambda_pre_raw[ : , 1] ~ normal(0, 3); // starting prior
  to_array_1d(logLambda_pre_raw[ : , 2 : maxDelay]) ~ normal(0, 1); // non-centered
  to_array_1d(logLambda_raw) ~ normal(0, 1); // non-centered
  
  // Delay dist
  //coefs for logit @ delay 0,..,maxDelay-1
  gamma ~ normal(mu_gamma, sd_gamma);
  // Changepoint effects
  eta_cp ~ normal(eta_cp_mu, eta_cp_sd);
  eta_extra ~ normal(eta_extra_mu, eta_extra_sd);
  
  // Prior for overdispersion
  //xi_negbinom ~ normal(0., 1.);
  
  // Random walk prior for percentage of cases with onset date
  sd_logitAlpha[ : ] ~ normal(0, 0.5); // half-normal due to constraint
  logitAlpha_raw[ : , 1] ~ normal(0, 2); // starting prior
  to_array_1d(logitAlpha_raw[ : , 2 : T]) ~ normal(0, 1); // non-centered
  
  // Likelihood
  for (j in 1 : J) {
    //for each case category j
    array[T] real expected_a = rep_array(0.1, T);
    // modeling phase (before time 1, no observations yet)
    for (t_pre in 1 : maxDelay) {
      for (d in (maxDelay - t_pre + 1) : maxDelay) {
        real expected_overall = exp(logLambda_pre[j, t_pre])
                                * p_pre[j, t_pre, d + 1];
        // add expected number of cases without observed onset arising from day t
        expected_a[t_pre + d - maxDelay] += expected_overall
                                            * (1 - alpha[j, 1]); //assume the same alpha value for all days in the modeling phase
      }
    }
    // inference (observations at time 1 to T)
    for (t in 1 : T) {
      for (d in 0 : min(T - t, maxDelay)) {
        real expected_overall = exp(logLambda[j, t]) * p[j, t, d + 1];
        // Likelihood for cases with observed onset
        rT[j, t, d + 1] ~ poisson(expected_overall * alpha[j, t] + 0.1); //, phi_negbinom);
        // add expected number of cases without observed onset arising from day t
        expected_a[t + d] += expected_overall * (1 - alpha[j, t]);
      }
    }
    // Likelihood for cases without observed onset
    aT[j,  : ] ~ poisson(expected_a); //, phi_negbinom);
  }
}
generated quantities {
  array[J, T - predLag, maxDelay + 1] int nT;
  array[J, T - predLag, maxDelay + 1] int n_aT;
  array[J, T - predLag] real ntInf;
  array[J, T - predLag] real n_atInf;
  array[J, T - predLag] real n_totalInf;
  
  for (j in 1 : J) {
    //for each case category j
    for (t in 1 : (T - predLag)) {
      for (d in 0 : maxDelay) {
        if ((t + d) <= T) {
          nT[j, t, d + 1] = rT[j, t, d + 1];
          n_aT[j, t, d + 1] = poisson_rng(exp(logLambda[j, t])
                                          * p[j, t, d + 1]
                                          * (1 - alpha[j, t]) + 0.1);
        } else {
          if (exp(logLambda[j, t]) * p[j, t, d + 1] > 1073740000) {
            print("j = ", j, "t = ", t, ", d = ", d, ", alpha = ", alpha[j, t], ", logLambda = ", logLambda[j], ", gamma = ", gamma, ", eta_cp = ", eta_cp, ", eta_extra = ", eta_extra, ", p =", p[j, t,  : ]);
          }
          nT[j, t, d + 1] = poisson_rng(exp(logLambda[j, t]) * p[j, t, d + 1]
                                        * alpha[j, t] + 0.1); //, phi_negbinom);
          n_aT[j, t, d + 1] = poisson_rng(exp(logLambda[j, t])
                                          * p[j, t, d + 1]
                                          * (1 - alpha[j, t]) + 0.1); //, phi_negbinom);
        }
      }
      //ntInf[t-(T-maxDelay)] = sum(nT[t-(T-maxDelay),]) - (t-(T-maxDelay)) * 0.1;
      ntInf[j, t] = sum(nT[j, t,  : ]);
      n_atInf[j, t] = sum(n_aT[j, t,  : ]);
      n_totalInf[j, t] = ntInf[j, t] + n_atInf[j, t];
    }
  }
}