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
  * Damped Holt's method as an innovation state space model, non-centered parameterization
  */
  vector holt_damped_process_noncentered(real alpha, real beta, real phi, real l_start, real b_start, vector noise_standardized, real noise_sd){
    int n = num_elements(noise_standardized);
    vector[n] b;
    vector[n] l;
    vector[n] y;
    real epsilon;
    
    epsilon = noise_standardized[1] * noise_sd;
    b[1] = phi * b_start + beta * epsilon;
    l[1] = l_start + phi * b_start + alpha * epsilon;
    y[1] = l_start + phi * b_start + epsilon;
    
    for(t in 2:n){
      epsilon = noise_standardized[t] * noise_sd;
      b[t] = phi * b[t-1] + beta * epsilon;
      l[t] = l[t-1] + phi * b[t-1] + alpha * epsilon;
      y[t] = l[t-1] + phi * b[t-1] + epsilon;
    }
    
    return(y);
  }
  
  /**
  * Cumulative product for column vector
  */
  vector cumulative_prod_column(vector x) {
    return exp(cumulative_sum(log(x)));
  }
  
  /**
  * Cumulative product of 1-x for column vector
  */
  vector cumulative_prod1m_column(vector x) {
    return exp(cumulative_sum(log1m(x)));
  }
  
  /**
  * Cumulative product for row vector
  */
  row_vector cumulative_prod_row(row_vector x) {
    return exp(cumulative_sum(log(x)));
  }
  
  /**
  * Cumulative product of 1-x for row vector
  */
  row_vector cumulative_prod1m_row(row_vector x) {
    return exp(cumulative_sum(log1m(x)));
  }
  
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
  
  vector softplus(vector x, real k) {
    return(log1p_exp(k * x) / k);
  }