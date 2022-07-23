/** ---------------------------------------------------------------------------
Sampling functions / transformations approximating count data distributions
---------------------------------------------------------------------------- */
 
  real approx_poisson_log_lpdf(vector y, vector mean_log) {
    int n = num_elements(y);
    vector[n] sigma2 = log1p_exp(-mean_log);
    vector[n] mu = mean_log - sigma2/2;
    return normal_lpdf(y | mu, sqrt(sigma2));
  }
  
  real approx_poisson_lpdf(vector y, vector mean_log) {
    int n = num_elements(y);
    vector[n] sigma2 = log1p_exp(-mean_log);
    vector[n] mu = mean_log - sigma2/2;
    return lognormal_lpdf(y | mu, sqrt(sigma2));
  }
  
  vector approx_poisson_lnorm_noncentered(vector mean_log, vector noise_standardized) {
    int n = num_elements(mean_log);
    vector[n] sigma2 = log1p_exp(-mean_log);
    vector[n] mu = mean_log - sigma2/2;
    return exp(mu + sqrt(sigma2) .* noise_standardized);
  }
  
  real approx_poisson_lnorm_noncentered(real mean_log, real noise_standardized) {
    return approx_poisson_lnorm_noncentered(rep_vector(mean_log,1), rep_vector(noise_standardized,1))[1];
  }