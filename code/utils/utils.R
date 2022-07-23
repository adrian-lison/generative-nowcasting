
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

weighted.median <- function(x, w) {
  w <- w[order(x)]
  x <- x[order(x)]

  prob <- cumsum(w) / sum(w)
  ps <- which.min(abs(prob - .5))
  return(x[ps])
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
  for (t in 1:T_all) { # iterating over reference dates
    for (d in (max(1, D - t + 2)):min(D + 1, T_all - t + 1)) { # iterating over delays
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
    triangle_rep[max(t - D, 1):min(t, T_all - D), T_all - t + 1] <- triangle[(max(1, D - t + 2)):min(D + 1, T_all - t + 1), t]
  }
  return(triangle_rep)
}

## ---------------------------------------------------------------
##                Discretized delay distributions               -
## ---------------------------------------------------------------

get_incubation_dist <- function(gamma_mean = 5.3, gamma_sd = 3.2, maxInc = 14) {
  # Incubation period probabilities are returned in forward order,
  # i.e. a period of zero days first and the longest period last
  gamma_shape <- (gamma_mean / gamma_sd)^2
  gamma_rate <- gamma_mean / (gamma_sd^2)
  # longest period (combines all periods >= maxInc)
  longest <- (1 - pgamma(maxInc, shape = gamma_shape, rate = gamma_rate))
  probs <- c(
    ddgamma(0:(maxInc - 1), shape = gamma_shape, rate = gamma_rate), # all except longest (discrete)
    longest
  )
  return(probs)
}

get_generation_dist <- function(gamma_mean = 4.8, gamma_sd = 2.3, maxGen = 10) {
  # The generation time probabilities are returned in forward order,
  # i.e. a generation time of one day first and the longest generation time last
  gamma_shape <- (gamma_mean / gamma_sd)^2
  gamma_rate <- gamma_mean / (gamma_sd^2)
  shortest <- pgamma(2, shape = gamma_shape, rate = gamma_rate)
  # longest period (combines all periods >= maxGen)
  longest <- (1 - pgamma(maxGen, shape = gamma_shape, rate = gamma_rate))
  probs <- c(
    shortest,
    ddgamma(2:(maxGen - 1), shape = gamma_shape, rate = gamma_rate), # all other (discrete)
    longest
  )
  return(probs)
}

get_discrete_lognormal <- function(meanlog, sdlog, maxX) {
  longest <- (1 - plnorm(maxX, meanlog = meanlog, sdlog = sdlog))
  probs <- c(
    sapply(0:(maxX - 1), function(x) plnorm(x + 1, meanlog = meanlog, sdlog = sdlog) - plnorm(x, meanlog = meanlog, sdlog = sdlog)), # all except longest (discrete)
    longest
  )
  return(probs)
}


get_hazard_from_p <- function(p) {
  p_cumulative <- c(0, cumsum(p))
  hazard <- p[1:length(p)] / (1 - p_cumulative[1:length(p)])
  hazard[is.na(hazard)] <- 0
  return(hazard)
}

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

get_prior <- function(variable, ...) {
  prior_list <- list()
  prior_values <- list(...)
  for (v in names(prior_values)) {
    prior_list[[paste0(variable, "_prior_", v)]] <- prior_values[[v]]
  }
  return(prior_list)
}

holt_damped_process_noncentered <- function(alpha, beta_star, phi, l_start, b_start, increments) {
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

## ---------------------------------------------------------------
##                        Cmdstanr functions                     -
## ---------------------------------------------------------------

#' Remove profiling statements from a character vector representing stan code
#'
#' @param s Character vector representing stan code
#'
#' @return A `character` vector of the stan code without profiling statements
remove_profiling <- function(s) {
  while (grepl("profile\\(.+\\)\\s*\\{", s, perl = TRUE)) {
    s <- gsub(
      "profile\\(.+\\)\\s*\\{((?:[^{}]++|\\{(?1)\\})++)\\}", "\\1", s,
      perl = TRUE
    )
  }
  return(s)
}

#' Write copies of the .stan files of a Stan model and its #include files
#' with all profiling statements removed.
#'
#' @param stan_file The path to a .stan file containing a Stan program.
#'
#' @param include_paths Paths to directories where Stan should look for files
#' specified in #include directives in the Stan program.
#'
#' @param target_dir The path to a directory in which the manipulated .stan
#' files without profiling statements should be stored. To avoid overriding of
#' the original .stan files, this should be different from the directory of the
#' original model and the `include_paths`.
#'
#' @return A `list` containing the path to the .stan file without profiling
#' statements and the include_paths for the included .stan files without
#' profiling statements
write_stan_files_no_profile <- function(stan_file, include_paths = NULL,
                                        target_dir = tempdir()) {
  # remove profiling from main .stan file
  code_main_model <- paste(readLines(stan_file, warn = FALSE), collapse = "\n")
  code_main_model_no_profile <- remove_profiling(code_main_model)
  if (!dir.exists(target_dir)) {
    dir.create(target_dir, recursive = T)
  }
  main_model <- cmdstanr::write_stan_file(
    code_main_model_no_profile,
    dir = target_dir,
    basename = basename(stan_file)
  )

  # remove profiling from included .stan files
  include_paths_no_profile <- rep(NA, length(include_paths))
  for (i in length(include_paths)) {
    include_paths_no_profile[i] <- file.path(
      target_dir, paste0("include_", i), basename(include_paths[i])
    )
    include_files <- list.files(
      include_paths[i],
      pattern = "*.stan", recursive = TRUE
    )
    for (f in include_files) {
      include_paths_no_profile_fdir <- file.path(
        include_paths_no_profile[i], dirname(f)
      )
      code_include <- paste(
        readLines(file.path(include_paths[i], f), warn = FALSE),
        collapse = "\n"
      )
      code_include_paths_no_profile <- remove_profiling(code_include)
      if (!dir.exists(include_paths_no_profile_fdir)) {
        dir.create(include_paths_no_profile_fdir, recursive = T)
      }
      cmdstanr::write_stan_file(
        code_include_paths_no_profile,
        dir = include_paths_no_profile_fdir,
        basename = basename(f)
      )
    }
  }
  return(list(model = main_model, include_paths = include_paths_no_profile))
}
