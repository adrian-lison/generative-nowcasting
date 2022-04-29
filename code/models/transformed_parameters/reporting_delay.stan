/** ---------------------------------------------------------------------------
Reporting delay model
---------------------------------------------------------------------------- */
  // temp variables for hazard conputation (order of dimensions reversed to avoid transposition)
    real contrib_occurrence_covariates;
    vector[n_delays] logit_haz_interval;
    vector[D] hazard;
    vector[D] cum_converse_haz; 
    
    // modeling phase
    for(t_pre in 1:D) {
      //delay distribution
      logit_haz_interval[1:n_delays] = gamma[1:n_delays] + W[t_pre,1:n_delays] * eta;
      hazard[1:D] = inv_logit(logit_haz_interval[delay_idx[1:D]]);
      
      cum_converse_haz = cumulative_prod_column(1-hazard[1:D]);
      p_pre[:, t_pre] = compute_prob_from_hazard(hazard,cum_converse_haz);
    }
    
    // inference phase
    for(t in 1:T) {
      //delay distribution
      contrib_occurrence_covariates = Z[D+t] * beta;
      logit_haz_interval[1:n_delays] = gamma[1:n_delays] + rep_vector(contrib_occurrence_covariates,n_delays) + W[D+t,1:n_delays] * eta;
      hazard[1:D] = inv_logit(logit_haz_interval[delay_idx[1:D]]);
      
      cum_converse_haz = cumulative_prod_column(1-hazard[1:D]);
      p[:, t] = compute_prob_from_hazard(hazard,cum_converse_haz);
    }