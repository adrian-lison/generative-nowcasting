# Fully Bayesian nowcasting of transmission dynamics

The nowcasting models in this repository are based on the general framework described in [Günther et al. 2021](https://doi.org/10.1002/bimj.202000112), i.e. they use a discrete-time survival model of the delay distribution (as originally proposed by [an der Heiden and Höhle 2014](https://doi.org/10.1111/biom.12194)) and an autocorrelated time series smoothing prior on the expected number of eventually observed events (as originally proposed by [McGough et al. 2020](https://doi.org/10.1371/journal.pcbi.1007735)).

The models in this repository have been developed further to account for cases with missing data and to integrate an underlying renewal process of infections. This allows to conduct missing data imputation, nowcasting and inference of transmission dynamics jointly in one model. Some of these aspects will also be implemented in [this](https://github.com/epiforecasts/epinowcast) package.

Detailed model specification and preprint coming!
