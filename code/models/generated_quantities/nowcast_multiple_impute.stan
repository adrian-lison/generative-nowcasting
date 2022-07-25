/** ---------------------------------------------------------------------------
Generated quantities for nowcast
---------------------------------------------------------------------------- */
  int nowcast_known[T];

  {
    int draw = discrete_range_rng(1, n_imputations);
    
    for (t in 1:(T)) {
      int nT_known[D+1];
      
      // occurrence dates where all events have been observed (D reached)
      nT_known[1:(1+min(D,T-t))] = reported_imputed[draw, t, 1:(1+min(D,T-t))];
      
      // occurrence dates where not all events have been observed yet
      if (D>T-t){
        if (overdispersion) {
          nT_known[(2+T-t):(D+1)] = neg_binomial_2_log_rng(lambda_log[n_lambda_pre+t]+p_log[(2+T-t):(D+1),t], phi_negbinom);
        } else {
          nT_known[(2+T-t):(D+1)] = poisson_log_rng(lambda_log[n_lambda_pre+t]+p_log[(2+T-t):(D+1),t]);
        }
      }
      
      // sum all events to the respective occurrence date
      nowcast_known[t] = sum(nT_known);
    }
  }