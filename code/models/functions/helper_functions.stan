/** ---------------------------------------------------------------------------
Helper functions for primitive operations
---------------------------------------------------------------------------- */

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
  * Softplus activation function
  */
  vector softplus(vector x, real k) {
    return(log1p_exp(k * x) / k);
  }
  
  /**
  * Element-wise log_sum_exp between two vectors
  */
  vector log_sum_exp_elementwise(vector x, vector y) {
    int n = num_elements(x); // we assume x and y have equal length
    vector[n] res;
    for (i in 1:n) {
      res[i] = log_sum_exp(x[i], y[i]);
    }
    return(res);
  }
  
  /**
  * Efficient dot product on log scale
  */
  real log_dot_product(vector x, vector y) {
    return(log_sum_exp(x + y));
  }
