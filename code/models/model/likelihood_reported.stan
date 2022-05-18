/** ---------------------------------------------------------------------------
Likelihood for reported events with known and unknown occurrence date
---------------------------------------------------------------------------- */
  // expected events with unknown occurrence date, by reporting date
  vector[T] expected_a=rep_vector(0.1,T);
  
  // modeling phase (before time 1, no observations yet)
  for (t_pre in 1:D){
    vector[t_pre] expected_overall = lambda[t_pre]*p_pre[(D-t_pre+2):(D+1),t_pre];
    // add expected number of events with unknown occurrence date arising from day t
    // assume the same alpha value for all days in the modeling phase
    expected_a[1:t_pre] += expected_overall[1:t_pre]*(1-alpha[1]);
  }
  
  // inference (observations at time 1 to T)
  for(t in 1:T) {
    vector[min(T-t,D)+1] expected_overall = lambda[D+t]*p[1:(min(T-t,D)+1),t];
    // likelihood for events with known occurrence date
    reported_known[t,1:(min(T-t,D)+1)] ~ neg_binomial_2(expected_overall*alpha[t]+0.1, phi_negbinom);
    // add expected number of events with unknown occurrence date arising from day t
    expected_a[t:(t+min(T-t,D))] += expected_overall*(1-alpha[t]);
  }
  
  // likelihood for events with unknown occurrence date
  reported_unknown ~ neg_binomial_2(expected_a, phi_negbinom);