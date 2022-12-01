/** ---------------------------------------------------------------------------
Reporting delay model
---------------------------------------------------------------------------- */
  // temp variables for hazard computation
    real contrib_occurrence_covariates;
    vector[n_delays] logit_haz_interval;
    vector[D] hazard;
    
    // modeling phase
    for(t_pre in 1:n_lambda_pre) {
      //delay distribution
      logit_haz_interval[1:n_delays] = gamma[1:n_delays] + W[t_pre,1:n_delays] * eta;
      hazard[1:D] = inv_logit(logit_haz_interval)[delay_idx[1:D]];
      
      p_log_pre[:, t_pre] = compute_prob_log_from_hazard(hazard[1:D]);
    }
    
    // inference phase
    for(t in 1:T) {
      //delay distribution
      contrib_occurrence_covariates = Z[n_lambda_pre+t] * beta;
      logit_haz_interval[1:n_delays] = gamma[1:n_delays] + rep_vector(contrib_occurrence_covariates,n_delays) + W[n_lambda_pre+t,1:n_delays] * eta;
      hazard[1:D] = inv_logit(logit_haz_interval)[delay_idx[1:D]];

      p_log[:, t] = compute_prob_log_from_hazard(hazard[1:D]);
    }
