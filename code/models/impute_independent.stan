data {
  int T; // number of days
  int D; // maximum reporting delay
  
  array[T, D+1] int reported_known; // reporting triangle
  array[T] int reported_unknown;
  
  // forward delay distribution
  vector[D+1] emp_forward_delay_dist;
}

parameters {

}

model {

}

generated quantities {
  array[T-D, D+1] int reported_unknown_imputed = rep_array(0, T-D, D+1); // imputed reporting triangle
  for (t in (D+1):T) {
    if (reported_unknown[t] > 0) {
      array[D+1] int backward = multinomial_rng(emp_forward_delay_dist, reported_unknown[t]);
      for (d in 1:min(D+1,t-D)) {
        reported_unknown_imputed[t - d + 1 - D, d] += backward[d];
      }
    }
  }
}

