/** ---------------------------------------------------------------------------
Functions for converting between probabilities and hazards
---------------------------------------------------------------------------- */
 
 /**
  * Compute probabilities from hazards
  */
  vector compute_prob_from_hazard(vector hazard) {
    int n = num_elements(hazard);
    vector[n+1] prob;
    vector[n] cum_converse_hazard;
    cum_converse_hazard = cumulative_prod1m_column(hazard[1:n]);
    prob[1] = hazard[1];
    prob[2:n] = hazard[2:n] .* cum_converse_hazard[1:(n-1)];
    prob[n+1] = cum_converse_hazard[n];
    return(prob);
  }
  
  /**
  * Compute log probabilities from hazards
  */
  vector compute_prob_log_from_hazard(vector hazard) {
    int n = num_elements(hazard);
    vector[n+1] prob_log;
    vector[n] cum_converse_hazard_log;
    cum_converse_hazard_log = cumulative_sum(log1m(hazard[1:n]));
    prob_log[1] = log(hazard[1]);
    prob_log[2:n] = log(hazard[2:n]) + cum_converse_hazard_log[1:(n-1)];
    prob_log[n+1] = cum_converse_hazard_log[n];
    return(prob_log);
  }