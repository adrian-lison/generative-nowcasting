/** ---------------------------------------------------------------------------
Helper functions
---------------------------------------------------------------------------- */
  /**
  * Non-centered and vectorized AR(1) process
  * Note that if the first value of the process is already included in the
  * increments_standardized vector, then you can provide a start_value = 0
  */
  vector ar1_process_noncentered_vec(real start_value, vector increments_standardized, real increment_sd) {
    return start_value + cumulative_sum(increments_standardized) * increment_sd;
  }
  
  /**
  * Cumulative product for column vector
  */
  vector cumulative_prod_column(vector x) {
    return exp(cumulative_sum(log(x)));
  }
  
  /**
  * Cumulative product for row vector
  */
  row_vector cumulative_prod_row(row_vector x) {
    return exp(cumulative_sum(log(x)));
  }
  
  /**
  * Cumulative product for row vector
  */
  vector compute_prob_from_hazard(vector hazard, vector cum_converse_hazard) {
    int n = num_elements(hazard);
    vector[n+1] prob;
    prob[1] = hazard[1];
    prob[2:n] = hazard[2:n] .* cum_converse_hazard[1:(n-1)];
    prob[n+1] = cum_converse_hazard[n];
    return(prob);
  }