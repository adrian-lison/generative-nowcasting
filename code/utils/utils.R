## ----------------------------------------------------------------
##                        Helper functions                       -
## ----------------------------------------------------------------

#' Clip values of a vector between a lower and an upper bound
#'
#' @param vec A vector with values to be clipped
#' @param LB Lower bound
#' @param UB Upper bound
#'
#' @return A vector with clipped/fenced values
fence <- function(vec, LB = -Inf, UB = Inf) {
  pmax(LB, pmin(vec, UB))
}

#' Compare entries of two vectors with potentially missing values
#'
#' @param v1 First vector
#' @param v2 Second vector, should be of same length as first one
#'
#' @return A boolean vector with element-wise comparisons. Comparisons including
#' NAs are coded as FALSE.
compareNA <- function(v1, v2) {
  same <- (v1 == v2) | (is.na(v1) & is.na(v2))
  same[is.na(same)] <- FALSE
  return(same)
}

#' Compute a weighted median of values
#'
#' @param x A vector with values to be aggregated
#' @param w Weights for each entry
#'
#' @return Weighted median of the values in x
weighted.median <- function(x, w) {
  w <- w[order(x)]
  x <- x[order(x)]

  prob <- cumsum(w) / sum(w)
  ps <- which.min(abs(prob - .5))
  return(x[ps])
}

#' Scale values proportionally between a minimum and maximum
#'
#' @param x A vector with values to scale
#' @param target_min The lower bound / minimum
#' @param target_max The upper bound / maximum
#'
#' @return A vector with the scaled values
scale_between <- function(x, target_min, target_max) {
  target_min + (x - min(x)) * ((target_max - target_min) / (max(x) - min(x)))
}

#' Takes a character vector of ranges (e.g. 5-23) and makes a comma-separated 
#' list of all values in the range
range_to_comma_list <- function(range_vector) {
  sapply(range_vector, function(x) {
    paste0(as.integer(str_split(x, "-")[[1]][1]):
             as.integer(str_split(x, "-")[[1]][2]),
           collapse = ","
           )
  })
}

#' Takes a reporting triangle by date of reference, and reshapes it into a
#' reporting triangle by date of report, excluding left-truncated dates.
#'
#' @param triangle DxT `matrix` representing a reporting triangle, i.e. columns
#' represent the date of reference and rows represent the delay
#'
#' @param D maximum delay
#'
#' @param pre_diff Difference between maximum delay and previously modeled
#' latent events
#'
#' @return A Dx(T-D) `matrix` representing a reporting triangle by date
#' of report. Columns represent the date of report (since first day that is
#' not left-truncated) and rows represent the delay.
#'
reporting_triangle_by_report <- function(triangle, D) {
  T_all <- dim(triangle)[2]
  triangle_rep <- matrix(0, nrow = D + 1, ncol = T_all - D)
  for (t in 1:T_all) { # iterat over reference dates
    for (d in (max(1, D - t + 2)):min(D + 1, T_all - t + 1)) { # iterate delays
      triangle_rep[d, t + d - 1 - D] <- triangle[d, t]
    }
  }
  return(triangle_rep)
}

#' Takes a reporting triangle by date of reference, and reshapes it into a
#' padded reporting triangle by date of report, excluding left-truncated dates.
#' This is a version optimized for stan, i.e. matrices are only accessed in
#' column-major order and the loop over delays has been vectorized. Note that
#' the dimensions and orientation of the returned matrix are different from the
#' non-optimized version.
#'
#' @param triangle DxT `matrix` representing a reporting triangle, i.e. columns
#' represent the date of reference and rows represent the delay
#'
#' @param D maximum delay
#'
#' @return A (T-D)xT `matrix` representing a reporting triangle by date
#' of report. Rows represent the date of report (since first day that is
#' not left-truncated) and columns represent the delay, but ragged. That is,
#' for a given reporting date, the sequence of reports with delay 1:(D+1) has
#' some leading and trailing zeros in the matrix. To obtain the overall cases
#' for a given reference date, one can then pre-multiply this matrix with a
#' T-length row vector of ones.
#'
reporting_triangle_by_report_padded <- function(triangle, D) {
  T_all <- dim(triangle)[2]
  triangle_rep <- matrix(0, nrow = T_all - D, ncol = T_all)
  for (t in 1:T_all) { # iterating over reference dates
    triangle_rep[max(t - D, 1):min(t, T_all - D), T_all - t + 1] <- triangle[
      (max(1, D - t + 2)):min(D + 1, T_all - t + 1), t
    ]
  }
  return(triangle_rep)
}

## ---------------------------------------------------------------
##                Discretized delay distributions               -
## ---------------------------------------------------------------

#' Get shape of a Gamma distribution given its mean and sd
get_gamma_shape_alternative <- function(gamma_mean, gamma_sd) {
  gamma_shape <- (gamma_mean / gamma_sd)^2
  return(gamma_shape)
}

#' Get rate of a Gamma distribution given its mean and sd
get_gamma_rate_alternative <- function(gamma_mean, gamma_sd) {
  gamma_rate <- gamma_mean / (gamma_sd^2)
  return(gamma_rate)
}

#' Get scale of a Gamma distribution given its mean and sd
get_gamma_scale_alternative <- function(gamma_mean, gamma_sd) {
  return(1 / get_gamma_rate_alternative(gamma_mean, gamma_sd))
}

#' Get PMF of a discretized Gamma distribution.
#'
#' This function accepts different parameterizations to specify the Gamma
#' distribution
#'
#' @param gamma_shape Shape parameter of the Gamma distribution
#' @param gamma_rate Rate parameter of the Gamma distribution.
#' @param gamma_scale Scale parameter of the Gamma distribution. Can be
#'   specified instead of the rate. Only has an effect if the rate is not
#'   specified.
#' @param gamma_mean Alternative parameterization: Mean of the Gamma
#' @param gamma_sd Alternative parameterization: Standard deviation of the Gamma
#' @param maxX Right truncation point
#' @param include_zero Should the distribution explicitly cover X=0, or should
#'   X=1 include the probability mass for X=0 too?
#' @param print_params Should the shape and rate parameters be printed?
#'
#' @return PMF of the discretized Gamma distribution
get_discrete_gamma <- function(gamma_shape,
                               gamma_rate,
                               gamma_scale,
                               gamma_mean,
                               gamma_sd,
                               maxX,
                               include_zero = T,
                               print_params = F) {
  if (missing(gamma_shape)) {
    if (missing(gamma_mean) || missing(gamma_sd)) {
      stop("No valid combination of parameters supplied", call. = F)
    }
    gamma_shape <- get_gamma_shape_alternative(gamma_mean, gamma_sd)
  }
  if (missing(gamma_rate)) {
    if (missing(gamma_scale)) {
      if (missing(gamma_mean) || missing(gamma_sd)) {
        stop("No valid combination of parameters supplied", call. = F)
      }
      gamma_rate <- get_gamma_rate_alternative(gamma_mean, gamma_sd)
    } else {
      gamma_rate <- 1 / gamma_scale
    }
  }

  # shortest period (combines periods 0 and 1)
  shortest <- pgamma(2, shape = gamma_shape, rate = gamma_rate)
  # longest period (combines all periods >= maxX)
  longest <- (1 - pgamma(maxX, shape = gamma_shape, rate = gamma_rate))

  if (include_zero) {
    probs <- c(
      # all except longest (discrete)
      ddgamma(0:(maxX - 1), shape = gamma_shape, rate = gamma_rate),
      longest
    )
  } else {
    probs <- c(
      shortest,
      # all other (discrete)
      ddgamma(2:(maxX - 1), shape = gamma_shape, rate = gamma_rate),
      longest
    )
  }

  if (print_params) {
    print(paste("Shape =", gamma_shape, "| Rate =", gamma_rate))
  }

  return(probs)
}

#' Get PMF of a discretized lognormal distribution.
#'
#' This function accepts both log-scale and unit-scale parameters to specify the
#' lognormal distribution
#'
#' @param meanlog Mean of log
#' @param sdlog Standard deviation of log
#' @param unit_mean Alternative parameterization: unit scale mean
#' @param unit_sd Alternative parameterization: unit scale sd
#' @param maxX Right truncation point
#' @param include_zero Should the distribution explicitly cover X=0, or should
#' X=1 include the probability mass for X=0 too?
#' @param print_params Should the log-level parameters be printed?
#'
#' @return PMF of the discretized lognormal
get_discrete_lognormal <- function(meanlog, sdlog, unit_mean = NULL,
                                   unit_sd = NULL, maxX, include_zero = T,
                                   print_params = F) {
  if (!is.null(unit_mean) && !is.null(unit_sd)) {
    sigma2 <- log((unit_sd / unit_mean)^2 + 1)
    meanlog <- log(unit_mean) - sigma2 / 2
    sdlog <- sqrt(sigma2)
  }
  # shortest period (combines periods 0 and 1)
  shortest <- plnorm(2, meanlog = meanlog, sdlog = sdlog)
  # longest period (combines all periods >= maxX)
  longest <- (1 - plnorm(maxX, meanlog = meanlog, sdlog = sdlog))

  if (include_zero) {
    probs <- c(
      sapply(0:(maxX - 1), function(x) {
        plnorm(x + 1, meanlog = meanlog, sdlog = sdlog) -
          plnorm(x, meanlog = meanlog, sdlog = sdlog)
      }), # all except longest (discrete)
      longest
    )
  } else {
    probs <- c(
      shortest,
      sapply(2:(maxX - 1), function(x) {
        plnorm(x + 1, meanlog = meanlog, sdlog = sdlog) -
          plnorm(x, meanlog = meanlog, sdlog = sdlog)
      }), # all other (discrete)
      longest
    )
  }

  if (print_params) {
    print(paste("meanlog =", meanlog, "| sdlog =", sdlog))
  }

  return(probs)
}

#' Convert a probability mass function to a hazard function
#'
#' @param p A vector representing the PMF
#'
#' @return The corresponding hazard function. Same length as PMF vector (last
#' element should be hazard = 1)
get_hazard_from_p <- function(p) {
  p_cumulative <- c(0, cumsum(p))
  hazard <- p[1:length(p)] / (1 - p_cumulative[1:length(p)])
  hazard[is.na(hazard)] <- 0
  return(hazard)
}

#' Convert a hazard function to a probability mass function
#'
#' @param hazard A vector with all discrete hazards
#'
#' @return A vector representing the corresponding PMF. This is one element
#' longer than the hazard vector (final hazard assumed to be 1).
get_p_from_hazard <- function(hazard) {
  cum_converse_hazard <- cumprod(1 - hazard)
  n <- length(hazard)
  p <- c(
    hazard[1],
    hazard[2:n] * cum_converse_hazard[1:(n - 1)],
    cum_converse_hazard[n]
  )
  return(p)
}

## ---------------------------------------------------------------
##                        State space models                     -
## ---------------------------------------------------------------

#' Simulate a trajectory from an innovations state space model
#' (Holt's linear trend method with dampening)
#'
#' @param alpha Smoothing parameter for level
#' @param beta_star Smoothing parameter for trend
#' @param phi Dampening factor
#' @param l_start Starting level
#' @param b_start Starting trend
#' @param increments Realized random innovations (epsilon)
#'
#' @return A vector with observations representing the realized trajectory
holt_damped_process_noncentered <- function(alpha, beta_star, phi, l_start,
                                            b_start, increments) {
  increments <- c(0, increments)
  n <- length(increments)
  beta <- alpha * beta_star

  if (phi == 0) {
    sum_b <- rep(0, n)
  } else if (phi == 1) {
    b <- b_start + beta * cumsum(c(0, increments))[1:n]
    sum_b <- cumsum(b)
  } else {
    b <- rep(NA, n)
    b[1] <- b_start
    for (t in 2:n) {
      b[t] <- phi * b[t - 1] + beta * increments[t - 1]
    }
    sum_b <- phi * cumsum(b)
  }

  y <- l_start + alpha * cumsum(c(0, increments))[1:n] + sum_b + increments

  return(y)
}
