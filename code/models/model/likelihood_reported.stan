/** ---------------------------------------------------------------------------
Likelihood for reported events with known and unknown occurrence date
---------------------------------------------------------------------------- */
  // expected events with unknown occurrence date, by reporting date
  int pre_diff = max(0, D-n_lambda_pre);
  matrix[D+1, n_lambda_pre+T] expected_a_log = rep_matrix(0, D+1, n_lambda_pre+T);
  
  // modeling phase (before time 1, no observations yet)
  for (t_pre in 1:n_lambda_pre){
    vector[t_pre] expected_overall_log = lambda_log[t_pre]+p_log_pre[(D-t_pre+2):(D+1), t_pre];
    // add expected number of events with unknown occurrence date arising from day t
    // assume the same alpha value for all days in the modeling phase
    expected_a_log[(D-t_pre+2):(D+1), t_pre] = expected_overall_log + alpha1m_log[1];
  }
  
  // inference (observations at time 1 to T)
  for(t in 1:T) {
    vector[1+min(T-t,D)] expected_overall_log = lambda_log[n_lambda_pre+t]+p_log[1:(min(T-t,D)+1),t];
    // likelihood for events with known occurrence date
    if (overdispersion) {
      reported_known[t,1:(1+min(T-t,D))] ~ neg_binomial_2_log(expected_overall_log+alpha_log[t], phi_negbinom);
    } else {
      reported_known[t,1:(1+min(T-t,D))] ~ poisson_log(expected_overall_log+alpha_log[t]);
    }
    // add expected number of events with unknown occurrence date arising from day t
    expected_a_log[1:(min(T-t,D)+1), n_lambda_pre+t] = expected_overall_log + alpha1m_log[t];
  }
  
  // reshape reporting triangle by date of report, leaving out left-truncated dates
  matrix[D+1, n_lambda_pre+T-D] expected_a_rep_log
    = reporting_triangle_by_report(expected_a_log, D);
  // compute log expected cases with missing occurrence for each date of report
  vector[n_lambda_pre+T-D] expected_a_rep_vec_log;
  for (i in 1:n_lambda_pre+T-D) {
    expected_a_rep_vec_log[i] = log_sum_exp(col(expected_a_rep_log,i));
  }
  
  // likelihood for events with unknown occurrence date
  if (overdispersion) {
    reported_unknown[(1+pre_diff):T] ~ neg_binomial_2_log(expected_a_rep_vec_log, phi_negbinom);
  } else {
    reported_unknown[(1+pre_diff):T] ~ poisson_log(expected_a_rep_vec_log);
  }
  
