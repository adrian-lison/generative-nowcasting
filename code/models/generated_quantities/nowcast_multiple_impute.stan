/** ---------------------------------------------------------------------------
Generated quantities for nowcast
---------------------------------------------------------------------------- */
  array[T] int nowcast_all;
  
  vector[T] predicted_known_rep = rep_vector(0,T); // model fit for cases with known onset date
  // note that predicted_known_rep is downward-biased for t<D+1 due to left-truncation
  {
    int draw = discrete_range_rng(1, n_imputations);
    
    for (t in 1:(T)) {
      array[D+1] int nT_known;
      
      // occurrence dates where all events have been observed (D reached)
      vector[1+min(D,T-t)] exp_cases_log_obs = lambda_log[n_lambda_pre+t]+p_log[1:(1+min(D,T-t)),t];
      if (overdispersion) {
        nT_known[1:(1+min(D,T-t))] = neg_binomial_2_log_rng(exp_cases_log_obs, phi_negbinom);
      } else {
        nT_known[1:(1+min(D,T-t))] = poisson_log_rng(exp_cases_log_obs);
      }
      
      
      // occurrence dates where not all events have been observed yet
      if (D>T-t){
        vector[D-(T-t)] exp_cases_log = lambda_log[n_lambda_pre+t]+p_log[(2+T-t):(D+1),t];
        if (overdispersion) {
          nT_known[(2+T-t):(D+1)] = neg_binomial_2_log_rng(exp_cases_log, phi_negbinom);
        } else {
          nT_known[(2+T-t):(D+1)] = poisson_log_rng(exp_cases_log);
        }
      }
      
      // assign cases with missing occurence date by date of report
      predicted_known_rep[t:min(t+D,T)] += to_vector(nT_known[1:min(1+D,1+T-t)]);
      
      // sum all events to the respective occurrence date
      nT_known[1:(1+min(D,T-t))] = reported_imputed[draw, t, 1:(1+min(D,T-t))]; // replace prediction for fully observed events with observation / imputed data
      nowcast_all[t] = sum(nT_known);
    }
  }