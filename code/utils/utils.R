
###########################################################################
###########################################################################
###                                                                     ###
###                          UTILITY FUNCTIONS                          ###
###                                                                     ###
###########################################################################
###########################################################################

## ----------------------------------------------------------------
##                        Helper functions                       -
## ----------------------------------------------------------------

fence <- function(vec, LB = -Inf, UB = Inf) pmax(LB, pmin(vec, UB))

compareNA <- function(v1, v2) {
  same <- (v1 == v2) | (is.na(v1) & is.na(v2))
  same[is.na(same)] <- FALSE
  return(same)
}

## ---------------------------------------------------------------
##                Discretized delay distributions               -
## ---------------------------------------------------------------

get_incubation_dist <- function(gamma_mean = 5.3, gamma_sd = 3.2, maxInc = 14) {
  # Note that the incubation period probabilities are returned in reverse order,
  # i.e. the longest period first and a period of zero days last
  gamma_shape <- (gamma_mean / gamma_sd)^2
  gamma_rate <- gamma_mean / (gamma_sd^2)
  # longest period (combines all periods >= maxInc)
  longest <- (1 - pgamma(maxInc, shape = gamma_shape, rate = gamma_rate))
  probs <- c(
    longest,
    ddgamma((maxInc - 1):0, shape = gamma_shape, rate = gamma_rate) # all other (discrete)
  )
  return(probs)
}

get_generation_dist <- function(gamma_mean = 4.8, gamma_sd = 2.3, maxGen = 10) {
  # Note that the generation time probabilities are returned in reverse order,
  # i.e. the longest generation time first and a generation time of one day last
  gamma_shape <- (gamma_mean / gamma_sd)^2
  gamma_rate <- gamma_mean / (gamma_sd^2)
  # longest period (combines all periods >= maxGen)
  longest <- (1 - pgamma(maxGen, shape = gamma_shape, rate = gamma_rate))
  shortest <- pgamma(2, shape = gamma_shape, rate = gamma_rate)
  probs <- c(
    longest,
    ddgamma((maxGen - 1):2, shape = gamma_shape, rate = gamma_rate), # all other (discrete)
    shortest
  )
  return(probs)
}
