/** ---------------------------------------------------------------------------
Delay distribution between latent events and occurrence events
---------------------------------------------------------------------------- */
  // maximum delay
  int L;

  // delay distribution
  // reversed, such that the probability for a delay of zero comes last
  vector[L+1] latent_delay_dist;
