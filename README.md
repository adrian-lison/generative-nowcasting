# Fully generative nowcasting of transmission dynamics

The nowcasting models in this repository are based on the general framework described in [Günther et al. 2021](https://doi.org/10.1002/bimj.202000112), i.e. they use a discrete-time survival model of the delay distribution (as originally proposed by [an der Heiden and Höhle 2014](https://doi.org/10.1111/biom.12194)) and an autocorrelated time series smoothing prior on the expected number of eventually observed events (as originally proposed by [McGough et al. 2020](https://doi.org/10.1371/journal.pcbi.1007735)).

The models in this repository have been developed further to account for cases with missing reference date and to integrate an underlying renewal process of infections. This allows to conduct missing data imputation, nowcasting and inference of transmission dynamics jointly in one model. Detailed model specification and preprint coming!

Over time, the extensions developed here will (partly) also be implemented in the [epinowcast](https://github.com/epiforecasts/epinowcast) package. In the following, a list of the current differences is provided.

### Incomplete data
- Generative model component for cases with missing reference date (share of missing dates modeled by date of reference) --> almost implemented in [epinowcast](https://github.com/epiforecasts/epinowcast), further testing needed
- Option to extend the latent time series of cases by reference date (and potentially the underlying infections and renewal process) further into the past. This allows to also include observed counts of cases with missing reference date that were reported further in the past (all the way until t=1, the first date in the reporting triangle). Without this extension, only counts after t=D can be used in an unbiased way. Note that the current implementation fixes the reporting delay distribution for t<1 to the distribution estimated for t=1.

### Reporting delay modeling
- Non-parametric model for reporting hazard, as in [Günther et al. 2021](https://doi.org/10.1002/bimj.202000112). Has been extended to both date of reference and date of report effect (already implemented in [epinowcast](https://github.com/epiforecasts/epinowcast). In turn, no parametric model for reporting hazard.
- Random walk or segmented change point model for reporting day effects. Both are already covered by the design matrix setup in [epinowcast](https://github.com/epiforecasts/epinowcast).
- Flexible reporting delay resolution which allows to group consecutive reporting delays together, assuming that they have the same hazard (identical baseline, and averaged hazard effects). Allows to shrink design matrices and reduce computations. This is useful in cases with very long-tailed reporting delays, where one can model delays individually for the first 30 days or so, and then gradually increase the window size. This may also help to relax pressure for choosing small maximum delays. Note that this is not fully pervasive in the implementation yet, i.e. the daily delay probabilities are still instantiated after the hazard effect modeling. Probably room for improvement.

### Expectation model
- Integrated renewal model with latent infections to estimate reproduction numbers and guide the nowcast. Renewal model is stochastic, i.e. realized infections are sampled. Generation interval and incubation period distributions are assumed as known and fixed. Simple random walk used to seed infections.
- Starting with the reference date time series, everything is on log scale. The renewal process however is not on the log scale.
- A few hacks for initializing the renewal model parameters already in a good spot such that sampling is faster.

### Time series priors (for expectation model)
- No time series prior, iid sampling (not recommended unless you have A LOT of data).
- Moving average model (simple, no parameters estimated, flexible window size). Vectorized.
- Innovations state space model (exponential smoothing): simple exponential smoothing, holt's linear trend method, holt's linear trend method with dampening. All parameters can either be estimated or fixed. Vectorized. Optional non-centered parameterization.
- All time series models can be differenced to an arbitrary order (by experience, this kills sampling performance however).

### Multiple imputation
This is mainly intended for comparison with the fully generative approach.
- Model to impute missing reference dates by estimating the backward delay distribution with a similar hazard model as in the nowcasting models. Can be used to get posterior draws for missing reference dates.
- Nowcasting model with an approximation to multiple imputation via likelihood averaging. Slightly faster than actual multiple imputation as it saves warm-up etc., but less precise.

### Post-processing utilities
- Validation of many nowcasts conducted over time using evaluation measures for point nowcasts, probabilistic nowcasts (scoring rules) and consistency of nowcasts. Plotting of nowcast for a certain delay over time, animated gifs with nowcasting over time.
- Simulation of line lists from a renewal process, with time-varying delay distributions.