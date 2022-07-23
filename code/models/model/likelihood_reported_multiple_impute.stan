/** ---------------------------------------------------------------------------
Likelihood for reported events with known occurrence date, multiple imputation
---------------------------------------------------------------------------- */

  // inference (observations at time 1 to T)
  for(t in 1:T) {
    vector[1+min(T-t,D)] expected_overall_log = lambda_log[n_lambda_pre+t]+p_log[1:(min(T-t,D)+1),t];
    // likelihood for events with known occurrence date
    if (overdispersion) {
      for (draw in 1:n_imputations) {
        target += 1/n_imputations * neg_binomial_2_log_lupmf(reported_imputed[draw,t,1:(1+min(T-t,D))] | expected_overall_log, phi_negbinom);
      }
    } else {
      for (draw in 1:n_imputations) {
      target += 1/n_imputations * poisson_log_lumpf( reported_imputed[draw,t,1:(1+min(T-t,D))] | expected_overall_log);
      }
    }
  }
  
