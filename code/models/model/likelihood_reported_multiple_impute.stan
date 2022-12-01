/** ---------------------------------------------------------------------------
Likelihood for reported events with known occurrence date, multiple imputation
---------------------------------------------------------------------------- */

  // inference (observations at time 1 to T)
  for(t in 1:T) {
    int n_obs = 1+min(T-t,D);
    vector[n_obs] expected_overall_log = lambda_log[n_lambda_pre+t]+p_log[1:n_obs,t];
    vector[n_imputations*n_obs] expected_overall_log_multiple;
    array[n_imputations*n_obs] int obs_multiple;
    for (draw in 1:n_imputations) {
      expected_overall_log_multiple[(n_obs*(draw-1)+1):(n_obs*draw)] = expected_overall_log;
      obs_multiple[(n_obs*(draw-1)+1):(n_obs*draw)] = reported_imputed[draw,t,1:n_obs];
    }
    // likelihood for events with known occurrence date
    if (overdispersion) {
        target += 1.0/n_imputations * neg_binomial_2_log_lupmf(obs_multiple | expected_overall_log_multiple, phi_negbinom);
    } else {
        target += 1.0/n_imputations * poisson_log_lupmf(obs_multiple | expected_overall_log_multiple);
    }
  }
  
