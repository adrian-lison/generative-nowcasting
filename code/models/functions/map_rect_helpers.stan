/** ---------------------------------------------------------------------------
Helper functions for stan's map_rect
---------------------------------------------------------------------------- */
 
  vector compute_logp(vector phi, vector theta, data array[] real x_r, data array[] int x_i){
    // unpack x_i
    int grainsize = x_i[1];
    int last_grainstep = x_i[2];
    int map_rect_group_size = x_i[3];
    int D = x_i[4];
    int n_delays = x_i[5];
    int n_beta = x_i[6];
    int n_eta = x_i[7];
    int delay_idx[D] = x_i[(7+1):(7+D)];
    
    // unpack theta
    vector[n_beta] beta = phi[1:n_beta];
    vector[n_eta] eta = phi[(n_beta+1):(n_beta+n_eta)];
    vector[n_delays] gamma = phi[(n_beta+n_eta+1):(n_beta+n_eta+n_delays)];
    
    row_vector[n_beta] Z;
    matrix[n_delays, n_eta] W;
    vector[n_delays] logit_haz_interval;
    vector[D] hazard;
    vector[(D+1) * grainsize] p_log;
    
    // unpack x_r
    real x_r_sub[map_rect_group_size];
    int sub_from, sub_to;
    for (grainstep in 1:last_grainstep) {
      sub_from = 1 + (grainstep-1) * map_rect_group_size;
      sub_to = grainstep * map_rect_group_size;
      x_r_sub = x_r[sub_from:sub_to];
      Z = to_row_vector(x_r_sub[1:n_beta]);
      W = to_matrix(x_r_sub[(n_beta + 1):(n_beta + (n_delays*n_eta))], n_delays, n_eta);
      // compute
      logit_haz_interval = gamma + rep_vector(Z * beta, n_delays) + W * eta;
      hazard = inv_logit(logit_haz_interval[delay_idx]);
      sub_from = 1 + (grainstep-1) * (D+1);
      sub_to = grainstep * (D+1);
      p_log[sub_from:sub_to] = compute_prob_log_from_hazard(hazard);
    }
    return(p_log);
  }