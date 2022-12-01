/** ---------------------------------------------------------------------------
Generated quantities for nowcast
---------------------------------------------------------------------------- */
  array[T] int nowcast_all;

  {
    int draw = discrete_range_rng(1, n_imputations);
    
    for (t in 1:(T)) {
      array[D+1] int nT_known;
      
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
      nowcast_all[t] = sum(nT_known);
    }
  }