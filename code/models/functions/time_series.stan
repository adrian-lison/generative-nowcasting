/** ---------------------------------------------------------------------------
Time series process functions
---------------------------------------------------------------------------- */
 
  /**
  * Non-centered and vectorized AR(1) process
  * Note that if the first value of the process is already included in the
  * increments_standardized vector, then you can provide a start_value = 0
  */
  vector random_walk(vector start_values, vector increments, int diff_order);
  
  vector random_walk(vector start_values, vector increments, int diff_order) {
    if (diff_order == 0) {
       return start_values[1] + cumulative_sum(append_row(0, increments));
     } else {
      vector[diff_order] next_start = start_values[2:(diff_order+1)];
      int next_n = num_elements(increments) + diff_order - 1;
      vector[next_n] diffs = random_walk(next_start, increments, diff_order - 1);
      return start_values[1] + cumulative_sum(append_row(0, diffs));
    }
  }
  
  /**
  * Damped Holt's method as an innovation state space model, non-centered parameterization, vectorized implementation
  * Linear method with additive noise and without seasonality
  */
  vector holt_damped_process(vector start_values, real alpha, real beta_star, real phi, vector noise, int diff_order);
  
  vector holt_damped_process(vector start_values, real alpha, real beta_star, real phi, vector noise, int diff_order) {
    if (diff_order == 0) {
      int n = num_elements(noise) + 1; // start value + length(noise)
      vector[n] epsilons = append_row(0, noise);
      vector[n] sum_b;
      real beta = alpha * beta_star; // below, we actually always need the product of alpha and beta
      
      if (phi == 0) {
        // special case: no trend
        sum_b = rep_vector(0, n);
      } else if (phi == 1) {
        // special case: trend, no dampening
        sum_b = cumulative_sum(start_values[2] + beta * cumulative_sum(append_row(0, epsilons))[1:n]);
      } else {
        // general case: trend and dampening
        vector[n] b = rep_vector(0, n);
        b[1] = start_values[2];
        // the following loop could also be implemented via matrix multiplication
        // but performance improvements are unlikely
        for (t in 2:n) {
          b[t] = phi * b[t - 1] + beta * epsilons[t - 1];
        }
        sum_b = phi * cumulative_sum(b);
      }
      return(start_values[1] + alpha * cumulative_sum(append_row(0, epsilons))[1:n] + sum_b + epsilons);
    } else {
      vector[diff_order+1] next_start = start_values[2:(diff_order+2)];
      int next_n = num_elements(noise) + diff_order;
      vector[next_n] diffs = holt_damped_process(next_start, alpha, beta_star, phi, noise, diff_order - 1);
      return start_values[1] + cumulative_sum(append_row(0, diffs));
    }
  }