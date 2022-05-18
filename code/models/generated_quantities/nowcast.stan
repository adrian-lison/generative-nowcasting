/** ---------------------------------------------------------------------------
Generated quantities for nowcast
---------------------------------------------------------------------------- */
  real nowcast_known[T];
  real nowcast_unknown[T];
  
  for(t in 1:(T)) {
    int nT_known[D+1];
    int nT_unknown[D+1];
    
    // occurrence dates where all events have been observed (D reached)
    nT_known[1:(1+min(D,T-t))] = reported_known[t, 1:(1+min(D,T-t))];
    nT_unknown[1:(1+min(D,T-t))] = neg_binomial_2_rng(lambda[D+t]*p[1:(1+min(D,T-t)),t]*(1-alpha[t])+0.1, phi_negbinom);
    
    // occurrence dates where not all events have been observed yet
    if(D>T-t){
      nT_known[(2+T-t):(D+1)] = neg_binomial_2_rng(lambda[D+t]*p[(2+T-t):(D+1),t]*alpha[t]+0.1, phi_negbinom);
      nT_unknown[(2+T-t):(D+1)] = neg_binomial_2_rng(lambda[D+t]*p[(2+T-t):(D+1),t]*(1-alpha[t])+0.1, phi_negbinom);
    }
    
    // sum all events to the respective occurrence date
    nowcast_known[t] = sum(nT_known);
    nowcast_unknown[t] = sum(nT_unknown);
  }