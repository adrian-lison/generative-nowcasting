/** ---------------------------------------------------------------------------
Delay distribution between latent events and occurrence events
---------------------------------------------------------------------------- */
  // maximum delay
  int L;

  // delay distribution
  // forward, the probability for a delay of zero comes first
  vector[L+1] latent_delay_dist;
