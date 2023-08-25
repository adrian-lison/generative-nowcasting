#' Create a result index data frame containing the file path and additional
#' information of each fit
#'
#' @param file_pattern Regex pattern that applies to both result and output
#'   files
#' @param result_folder Folder in which nowcasting results are stored
#' @param result_subfolder Path to the folder which contains all output and
#'   result files
#' @param resultfile_pattern Regex pattern for the result files
#' @param outputfile_pattern Regex pattern for the output files
#' @param result_info Should additional info on nowcasts be collected (e.g.
#'   maximum delay and nowcast date)?
#'
#' @returns A `data.frame` with columns: id, now, maxDelay, result_file,
#'   output_file
define_results <- function(file_pattern, result_folder, result_subfolder,
                           resultfile_pattern = NULL, outputfile_pattern = NULL,
                           result_info = T) {
  # determine result file pattern
  if (is.null(resultfile_pattern)) {
    resultfile_pattern <- paste0(file_pattern, "_\\d+_result\\.rds")
  } else {
    resultfile_pattern <- paste0(file_pattern, resultfile_pattern)
  }
  if (is.null(outputfile_pattern)) {
    outputfile_pattern <- paste0(file_pattern, "_\\d+_.*\\.out")
  } else {
    outputfile_pattern <- paste0(file_pattern, outputfile_pattern)
  }
  if (missing(result_subfolder)) result_subfolder <- file_pattern

  # get file paths
  result_path <- file.path(result_folder, result_subfolder)
  files_result <- data.frame(result_file = file.path(result_path, list.files(
    path = here::here(result_path), pattern = resultfile_pattern, full.names = F
  )))
  files_output <- data.frame(output_file = file.path(result_path, list.files(
    path = here::here(result_path), pattern = outputfile_pattern, full.names = T
  )))

  # determine IDs
  resultfile_start_end <- str_split(resultfile_pattern, pattern = "\\\\d\\+")[[1]]
  extract_result_id_pattern <- paste0(
    "(?<=", resultfile_start_end[1], ")\\d+(?=", resultfile_start_end[2], ")"
  )
  files_result$id <- str_extract(
    files_result$result_file, extract_result_id_pattern
  )
  files_result$result_date <- as.Date(str_extract(
    files_result$result_file, "(?<=for_)\\d{4}-\\d{2}-\\d{2}"
  ))

  outputfile_start_end <- str_split(
    outputfile_pattern,
    pattern = "\\\\d\\+"
  )[[1]]

  extract_output_id_pattern <- paste0(
    "(?<=", outputfile_start_end[1],
    ")\\d+(?=", outputfile_start_end[2], ")"
  )

  files_output$id <- str_extract(
    files_output$output_file,
    extract_output_id_pattern
  )

  if (any(duplicated(files_result$id))) {
    warning("There are duplicated IDs in the result files.")
  }

  results <- files_result %>%
    full_join(files_output, by = "id")

  # get additional info on individual result files
  if (result_info) {
    additional_info <- lapply(results$result_file, function(x) {
      try(
        {
          result_f <- read_rds(here::here(x))
          if ("now" %in% names(result_f$summary)) {
            return(list(
              result_f$summary$now,
              result_f$summary$fit_datetime,
              result_f$summary$D,
              result_f$dataset_index,
              result_f$date_index
            ))
          } else {
            return(list(NA, NA, NA, NaN, NaN))
          }
        },
        silent = T
      )
      return(list(NA, NA, NA, NaN, NaN))
    })

    results$now <- as.Date(sapply(additional_info, function(x) x[[1]]))
    results$fit_datetime <- sapply(additional_info, function(x) x[[2]])
    results$maxDelay <- sapply(additional_info, function(x) x[[3]])
    results$dataset_index <- sapply(additional_info, function(x) x[[4]])
    results$date_index <- sapply(additional_info, function(x) x[[5]])
  } else {
    results$now <- NA
    results$fit_datetime <- NA
    results$maxDelay <- NA
    results$dataset_index <- NaN
    results$date_index <- NaN
  }

  if (all(is.list(results$dataset_index))) {
    results$dataset_index <- NaN
  }
  if (all(is.list(results$date_index))) {
    results$date_index <- NaN
  }

  results <- results %>%
    relocate(id, dataset_index, date_index, now, fit_datetime, maxDelay,
      .before = 1
    ) %>%
    mutate(across(c(id, dataset_index, date_index), as.integer)) %>%
    arrange(id)

  # remove duplicates
  if (any(duplicated(files_output$id))) {
    warning("There are duplicated IDs in the output files.")
  }
  results <- results %>%
    arrange(desc(fit_datetime), desc(output_file)) %>%
    dplyr::distinct(across(-c(output_file))) # take the latest version

  return(results)
}

#' Get result indices for a list of nowcasting results from different models
#' stored in separate folders.
#'
#' @param result_folder Folder in which nowcasting results are stored
#' @param result_info Should additional info on nowcasts be collected (e.g.
#'   maximum delay and nowcast date)?
#' @param list_of_sources Named list with different result subfolders as values.
#'   Can be used to collect results from different nowcasting job arrays, e.g.
#'   different nowcasting models fitted on the same data.
#' @param overwrite Should the result indices stored in the corresponding result
#'   subfolders be overwritten or reused? Should be set to TRUE whenever results
#'   were updated.
#' @param resultfile_pattern Regex pattern for the result files
#' @param outputfile_pattern Regex pattern for the output files
#'
#' @returns A named `list` corresponding to `list_of_sources`, with the
#'   different result indices.
define_result_list <- function(result_folder = "results", result_info = T, list_of_sources,
                               overwrite = F, resultfile_pattern = NULL,
                               outputfile_pattern = NULL) {
  res_list <- as.list(list_of_sources)

  for (i in seq_along(res_list)) {
    res_file <- here::here(result_folder, res_list[[i]], "results_index.rds")
    if (!overwrite && file.exists(res_file)) {
      res_list[[i]] <- readRDS(res_file)
    } else {
      print(paste("Computing result index for", names(res_list)[i]))
      res_list[[i]] <- define_results(
        file_pattern = res_list[[i]],
        result_folder = result_folder,
        result_info = result_info,
        resultfile_pattern = resultfile_pattern,
        outputfile_pattern = outputfile_pattern
      )
      saveRDS(res_list[[i]], res_file)
    }
  }

  maxDelay <- max(sapply(res_list, function(x) max(x$maxDelay, na.rm = T)))
  if (any(sapply(res_list, function(x) any(x$maxDelay != maxDelay, na.rm = T)))) {
    warning("Not all nowcasts have the same maximum delay.")
  }

  return(res_list)
}

#' Load nowcasting results from simulated data alongside error metrics
#'
#' @param res_list Result indices for different models, see
#'   `define_result_list()`.
#' @param maxDelay Maximum assumed delay used by the models.
#' @param reference_date The starting reference date to map the simulation time
#'   points to dates.
#' @param ground_truth_sim_list A `list` of simulated ground truth values, each
#'   entry in the list represents a different simulation run.
#' @param keep_posterior Should the posterior samples from the fitted models
#'   also be returned or only used to compute performance metrics and then
#'   dropped?
#' @param include_naive_nowcasts Should naive nowcasts (no adjustment for right
#'   truncation) of cases be included?
#' @param overwrite Should the results stored in the corresponding result
#'   subfolders be overwritten or reused? Should be set to TRUE whenever results
#'   were updated.
#' @param mc_cores The loading of results can be parallelized to several cores.
#'
#' @return A named `list` corresponding to `res_list`, with the different
#'   results and error metrics.
load_results_sim <- function(res_list, maxDelay, reference_date,
                             ground_truth_sim_list, keep_posterior = TRUE,
                             include_naive_nowcasts = TRUE,
                             overwrite = FALSE,
                             mc_cores = 4) {
  results_sim <- lapply(seq_along(res_list), function(i) {
    parent_dir <- dirname(res_list[[i]]$result_file[
      min(which(!is.na(res_list[[i]]$result_file)))
    ])
    res_file <- here::here(parent_dir, "all_results.rds")

    if (!overwrite && file.exists(res_file)) {
      print(paste("Loading already existing results for", names(res_list)[i]))
      res <- readRDS(res_file)
    } else {
      print(paste("Extracting results for", names(res_list)[i]))

      res <- get_nowcasts(
        res_list[[i]],
        maxDelay,
        reference_date,
        ground_truth_sim = ground_truth_sim_list,
        keep_posterior = keep_posterior,
        mc_cores = mc_cores
      )

      if (include_naive_nowcasts) {
        print("--> Including naive nowcasts")
        res <- add_naive_nowcasts_simulated(
          res, ground_truth_sim_list, reference_date
        )
      }

      print("--> Saving")
      saveRDS(res, res_file)

      return(res)
    }
  })

  names(results_sim) <- names(res_list)
  return(results_sim)
}

#' Load nowcasting results from simulated data alongside error metrics
#'
#' @param res_list Result indices for different models, see
#'   `define_result_list()`.
#' @param maxDelay Maximum assumed delay used by the models.
#' @param reference_date The first observed reference date.
#' @param line_list_empirical If available, observed or approximated ground
#'   truth values.
#' @param keep_posterior Should the posterior samples from the fitted models
#'   also be returned or only used to compute performance metrics and then
#'   dropped?
#' @param include_naive_nowcasts Should naive nowcasts (no adjustment for right
#'   truncation) of cases be included?
#' @param overwrite Should the results stored in the corresponding result
#'   subfolders be overwritten or reused? Should be set to TRUE whenever results
#'   were updated.
#' @param mc_cores The loading of results can be parallelized to several cores.
#' @param consolidation_lags At what lags should nowcasts be treated as
#'   consolidated? These will be used to estimate approximate ground truth
#'   values for the empirical data.
#'
#' @return A named `list` corresponding to `res_list`, with the different
#'   results and error metrics.
load_results_real <- function(res_list, maxDelay, reference_date,
                              line_list_empirical, keep_posterior = FALSE,
                              include_naive_nowcasts = TRUE, 
                              overwrite = FALSE,
                              mc_cores = 1, consolidation_lags = 7:14) {
  # Run once without ground truth emp
  n_res_list <- lapply(res_list, function(res) {
    gc()
    all_results_path <- here::here(dirname(res$result_file[1]), "all_results.rds")
    if (!overwrite & file.exists(all_results_path)) {
      print(paste(all_results_path, "already exists"))
      return(readRDS(all_results_path))
    } else {
      ncres <- get_nowcasts(res, maxDelay, reference_date,
        keep_posterior = F,
        mc_cores = mc_cores
      )
      return(ncres)
    }
  })

  if (!is.null(line_list_empirical)) {
    # Get ground empirical ground truth (proxy)
    ground_truth_emp <- get_ground_truth_empirical(
      bind_rows(n_res_list, .id = "model") %>%
        filter(R_model == "renewal", model %in% c("fully_generative")),
      line_list_empirical, maxDelay,
      consolidation_lags = consolidation_lags
    )
    
    # Run again with ground empirical truth and save
    n_res_list <- lapply(res_list, function(res) {
      gc()
      all_results_path <- here::here(dirname(res$result_file[1]), "all_results.rds")
      if (!overwrite & file.exists(all_results_path)) {
        res_load <- readRDS(all_results_path)
        if (any(grepl("_true", names(res_load)))) {
          print(paste(all_results_path, "already exists with ground truth"))
          return(res_load)
        }
      }
      ncres <- get_nowcasts(
        res,
        maxDelay,
        reference_date,
        keep_posterior = keep_posterior,
        mc_cores = mc_cores,
        ground_truth_emp = ground_truth_emp
      )
      saveRDS(ncres, all_results_path)
      return(ncres)
    })
    
    # Add naive nowcasts
    if (include_naive_nowcasts) {
      naive_nc <- get_naive_nowcasts_empirical(
        line_list_empirical, as.integer(max(n_res_list[[1]]$delay, na.rm = T))
      )
      n_res_list <- lapply(n_res_list, function(x) {
        x %>% left_join(naive_nc, by = c("nowcast_date", "date" = "event1_date"))
      })
    }
  }
  
  return(n_res_list)
}

#' Get diagnostic information of each fit, along with some summaries
#'
#' @return A `list` with summary data frames and the full diagnostic data frames
get_diagnostics <- function(results, n_chains = 4, overwrite = FALSE, horizon = 4*7) {
  diags_file <- here::here(dirname(results$result_file[
    min(which(!is.na(results$result_file)))
  ]), "results_diagnostics.rds")

  if (!overwrite && file.exists(diags_file)) {
    print("Loading preprocessed diagnostics.")
    diags <- readRDS(diags_file)
    return(diags)
  } else {
    print("Extracting diagnostics.")

    diags <- bind_rows(apply(results, 1, function(res) {
      if (is.na(res["result_file"])) {
        return(data.frame(id = res["id"]))
      }
      nowcast <- read_rds(here::here(res["result_file"]))
      if (!("diagnostic_summary" %in% names(nowcast))) {
        return(data.frame(id = res["id"]))
      }

      diags <- nowcast$diagnostic_summary %>% as.data.frame()
      diags$id <- res["id"]
      diags$fit_datetime <- nowcast$summary$fit_datetime
      diags$dataset_index <- nowcast$dataset_index
      diags$date_index <- nowcast$date_index
      diags$now <- nowcast$summary$now
      
      ess_params <- str_trim(str_extract(nowcast$diagnostics$stdout,
        pattern = "(?<=effective draws per transition:\\n).*(?=\\n)"
      ))
      ess_nowcast <- as.integer(str_extract_all(
        ess_params, "(?<=nowcast_all\\[)\\d+(?=\\])")[[1]])
      diags$ess_nowcast <- list(ess_nowcast[
        with(nowcast$summary,ess_nowcast > T+L-horizon)
        ])
      ess_R <- as.integer(str_extract_all(
        ess_params, "(?<=R\\[)\\d+(?=\\])")[[1]])
      diags$ess_R <- list(ess_R[
        with(nowcast$summary,ess_R > T+L-max_gen-horizon)
        ])
      
      rhat_params <- str_trim(str_extract(nowcast$diagnostics$stdout,
        pattern = "(?<=R-hat greater than \\d[[:punct:]]\\d\\d:\\n).*(?=\\n)"
      ))
      rhat_nowcast <- as.integer(str_extract_all(
        rhat_params, "(?<=nowcast_all\\[)\\d+(?=\\])")[[1]])
      diags$rhat_nowcast <- list(rhat_nowcast[
        with(nowcast$summary,rhat_nowcast > T+L-horizon)
      ])
      rhat_R <- as.integer(str_extract_all(
        rhat_params, "(?<=R\\[)\\d+(?=\\])")[[1]])
      diags$rhat_R <- list(rhat_R[
      with(nowcast$summary,rhat_R > T+L-max_gen-horizon)
      ])
      return(diags)
    }))

    diags <- dplyr::relocate(diags, c(
      id, fit_datetime, dataset_index, date_index, now
      )) %>%
      group_by(id, fit_datetime, dataset_index, date_index, now) %>%
      summarize(
        chains_divergent = sum(num_divergent > 40),
        chains_max_treedepth = sum(num_max_treedepth > 0),
        chains_low_ebfmi = sum(ebfmi < 0.3),
        across(c(ess_nowcast, ess_R, rhat_nowcast, rhat_R), 
               function(x) paste(unique(unlist(x)), collapse = ",")),
        .groups = "drop"
      ) %>% 
      mutate(across(c(ess_nowcast, ess_R, rhat_nowcast, rhat_R),
                    function(x) ifelse(x == "NA", "", x)))

    saveRDS(diags, diags_file)
    return(diags)
  }
}

#' Get information on how many cases longer than the maxDelay where truncated 
#' in each fit
get_truncation_info <- function(results) {
  truncs <- bind_rows(apply(results, 1, function(res) {
    if (is.na(res["output_file"])) {
      return(data.frame(id = res["id"]))
    }
    return(data.frame(
      id = res["id"],
      truncated = str_extract(read_file(res["output_file"]),
        pattern = "(?<=Removed )\\d+(?= cases with delay > \\d+\\.)"
      )
    ))
  })) %>%
    mutate(truncated = as.integer(truncated))
  rownames(truncs) <- NULL

  n_fits_truncated <- nrow(truncs %>% filter(truncated > 0))

  summary_truncation <- summary(truncs %>%
    filter(truncated > 0) %>%
    pull(truncated))

  return(list(
    n_fits_truncated = n_fits_truncated,
    summary_truncation = summary_truncation,
    truncs = truncs
  ))
}

#' Get a vector of dates for which a valid fit is missing
get_missing_dates <- function(results, all_dates = NULL) {
  if (is.null(all_dates)) {
    all_dates <- seq.Date(
      min(results$now, na.rm = T),
      max(results$now, na.rm = T),
      by = "1 day"
    )
  }
  missing_dates <- results %>%
    group_by(dataset_index) %>%
    summarize(
      n_missing = length(setdiff(all_dates, now)),
      missing_dates_list = list(sort(as.Date(setdiff(all_dates, now)))),
      missing_dates = paste(sapply(missing_dates_list, as.Date), collapse = ", "),
      .groups = "keep"
    )
  return(missing_dates)
}

#' Load all nowcasting results from files and store them in a common data.frame
#' @param results Result index for a certain model, see
#'   `define_results()`.
#' @param maxDelay Maximum assumed delay used by the model.
#' @param reference_date The first observed or simulated reference date.
#' @param ground_truth_sim A `data.frame` with simulated ground truth values.
#' @param ground_truth_emp A `data.frame` with observed or approximated ground
#' truth values.
#' @param keep_posterior Should the posterior samples from the fitted models
#'   also be returned or only used to compute performance metrics and then
#'   dropped?
#' @param mc_cores The loading of results can be parallelized to several cores.
get_nowcasts <- function(results, maxDelay, reference_date = NULL,
                         ground_truth_sim = NULL, ground_truth_emp = NULL,
                         keep_posterior = TRUE, mc_cores = 4) {
  all_res <- parallel::mclapply(
    apply(results, 1, function(res) res, simplify = F), function(res) {
      gc()
      verbose <- F

      if (is.na(res["result_file"]) | is.na(res["now"])) {
        return(data.frame(
          id = res["id"],
          nowcast_date = as.Date(res["now"]),
          dataset_index = res["dataset_index"],
          date_index = res["date_index"]
        ))
      }

      nowcast_file <- read_rds(here::here(res["result_file"]))
      if (!("summary" %in% names(nowcast_file))) {
        return(data.frame(
          id = res["id"],
          nowcast_date = as.Date(res["now"]),
          dataset_index = res["dataset_index"],
          date_index = res["date_index"]
        ))
      }

      if (verbose) print("Nowcast")
      if (!(all(c(
        "nowcast_all",
        "predicted_all_rep"
      ) %in% names(nowcast_file$summary$nowcast)))) {
        if (any(c(
          "nowcast",
          "nowcast_known",
          "nowcast_unknown",
          "predicted_known_rep"
        ) %in% names(nowcast_file$summary$nowcast))) {
          nowcast_file$summary$nowcast <- add_nowcast_all(nowcast_file$summary)
        }
      }

      if ("nowcast" %in% names(nowcast_file$summary)) {
        # read Nt nowcasts
        nc <- nowcast_file$summary$nowcast %>%
          filter(.width == 0.95) %>%
          select(-.width)
      } else {
        nc <- data.frame() %>% transmute(date = as.Date(c()))
      }

      if (verbose) print("R")
      # read Rt nowcasts
      nc_R_list <- list()
      if ("R" %in% names(nowcast_file$summary)) {
        nc_R_list[["default"]] <- nc %>% full_join(
          nowcast_file$summary$R %>%
            rename(R.lower = .lower, R.upper = .upper) %>%
            mutate(date = as.Date(as.character(date))),
          by = "date"
        )
      }
      if ("R_epiestim" %in% names(nowcast_file$summary)) {
        nc_R_list[["epiestim"]] <- nc %>% full_join(
          nowcast_file$summary$R_epiestim %>%
            rename(R.lower = .lower, R.upper = .upper) %>%
            mutate(date = as.Date(as.character(date))),
          by = "date"
        )
      }
      if ("R_renewal" %in% names(nowcast_file$summary)) {
        nc_R_list[["renewal"]] <- nc %>% full_join(
          nowcast_file$summary$R_renewal %>%
            rename(R.lower = .lower, R.upper = .upper) %>%
            mutate(date = as.Date(as.character(date))),
          by = "date"
        )
      }
      nc <- bind_rows(nc_R_list, .id = "R_model")
      remove(nc_R_list)

      if (verbose) print("delay")
      # read mean estimated delays
      if ("mean_delay" %in% names(nowcast_file$summary)) {
        nc <- nc %>% full_join(
          nowcast_file$summary$mean_delay %>%
            rename(mean_delay.lower = .lower, mean_delay.upper = .upper),
          by = "date"
        )
      }

      if (verbose) print("missingness")
      # read estimated proportion of known reference dates
      if ("fraction_complete" %in% names(nowcast_file$summary)) {
        nc <- nc %>% full_join(
          nowcast_file$summary$fraction_complete %>%
            rename(alpha.lower = .lower, alpha.upper = .upper),
          by = "date"
        )
      }

      if (verbose) print("posterior nowcast")
      # read posterior samples of nowcast
      if ("posterior_nowcast" %in% names(nowcast_file$summary)) {
        if (all(c("nowcast_known", "nowcast_unknown") %in% 
                names(nowcast_file$summary$posterior_nowcast)) &
          !("nowcast_all" %in% names(nowcast_file$summary$posterior_nowcast))) {
          nowcast_file$summary$posterior_nowcast <-
            nowcast_file$summary$posterior_nowcast %>%
            mutate(nowcast_all = nowcast_known + nowcast_unknown, 
                   .after = "nowcast_unknown")
        }

        nc <- nc %>%
          left_join(
            nowcast_file$summary$posterior_nowcast %>%
              group_by(date) %>%
              summarize(across(any_of(
                c(
                  "nowcast_known",
                  "nowcast_unknown",
                  "nowcast_all",
                  "predicted_missing_rep",
                  "predicted_known_rep"
                )
              ), list)) %>%
              plyr::rename(
                c(
                  "nowcast_known" = "nowcast_known_posterior",
                  "nowcast_unknown" = "nowcast_unknown_posterior",
                  "nowcast_all" = "nowcast_all_posterior",
                  "predicted_missing_rep" = "predicted_missing_rep_posterior",
                  "predicted_known_rep" = "predicted_known_rep_posterior"
                ),
                warn_missing = F
              ),
            by = "date"
          )
      }

      if (verbose) print("posterior R")
      # read posterior samples of R nowcasts
      nc_R_list <- list()
      if ("posterior_R" %in% names(nowcast_file$summary)) {
        nc_R_list[["default"]] <- nowcast_file$summary$posterior_R %>%
          group_by(date) %>%
          summarize(across(R, list)) %>%
          rename(R_posterior = R) %>%
          mutate(date = as.Date(as.character(date)))
      }
      if ("posterior_R_epiestim" %in% names(nowcast_file$summary)) {
        nc_R_list[["epiestim"]] <- nowcast_file$summary$posterior_R_epiestim %>%
          group_by(date) %>%
          summarize(across(R, list)) %>%
          rename(R_posterior = R) %>%
          mutate(date = as.Date(as.character(date)))
      }
      if ("posterior_R_renewal" %in% names(nowcast_file$summary)) {
        nc_R_list[["renewal"]] <- nowcast_file$summary$posterior_R_renewal %>%
          group_by(date) %>%
          summarize(across(R, list)) %>%
          rename(R_posterior = R) %>%
          mutate(date = as.Date(as.character(date)))
      }

      if (length(nc_R_list) > 0) {
        nc <- nc %>% left_join(
          bind_rows(nc_R_list, .id = "R_model"),
          by = c("R_model", "date")
        )
      }
      remove(nc_R_list)

      if (verbose) print("identifiers")
      # add identifiers
      nc <- nc %>%
        mutate(
          id = res["id"],
          nowcast_date = res["now"],
          dataset_index = res["dataset_index"],
          date_index = res["date_index"],
          .before = date
        ) %>%
        mutate(
          nowcast_date = as.Date(nowcast_date),
          delay = nowcast_date - date,
          .before = date
        ) %>%
        mutate(across(where(is.character), trimws))

      if (verbose) print("ground truth")
      # add ground truth
      if (!(is.null(ground_truth_sim))) {
        nc <- nc %>% add_ground_truth_simulated(
          ground_truth_sim, maxDelay, reference_date
        )
      }
      if (!(is.null(ground_truth_emp))) {
        nc <- nc %>% left_join(ground_truth_emp, by = "date")
      }

      # add performance metrics
      if (verbose) print("performance")
      nc <- nc %>% add_performance()

      # remove posterior if specified
      if (!keep_posterior) {
        nc <- nc %>% select(-contains("_posterior"))
      }
      return(nc)
    },
    mc.cores = mc_cores
  )

  gc()

  print("--> Merging all results.")
  all_res <- bind_rows(all_res)

  return(all_res)
}

# If missing, add nowcasts for all cases, i.e. the sum of cases with known and
# unknown reference date, (same for predicted cases) in hindsight (assuming that
# posterior nowcast is available). This returns the nowcast$summary$nowcast
# data.frame
add_nowcast_all <- function(nowcast_summary) {
  if (!("posterior_nowcast" %in% names(nowcast_summary))) {
    stop("Cannot add nowcast_all without posterior_nowcast.")
  }

  # nowcast all
  if (!("nowcast_all" %in% names(nowcast_summary$nowcast))) {
    if (all(c("nowcast_known", "nowcast_unknown") %in%
      names(nowcast_summary$nowcast))) {
      nowcast_summary$posterior_nowcast$nowcast_all <-
        nowcast_summary$posterior_nowcast$nowcast_known +
        nowcast_summary$posterior_nowcast$nowcast_unknown
    } else if (("nowcast_known" %in% names(nowcast_summary$nowcast)) &
      !("nowcast_unknown" %in% names(nowcast_summary$nowcast))) {
      # assume that there were not missing reference dates, hence "known" and "all" are identical
      nowcast_summary$posterior_nowcast$nowcast_all <-
        nowcast_summary$posterior_nowcast$nowcast_known
    }
  }

  # predicted all
  if (!("predicted_all_rep" %in% names(nowcast_summary$nowcast))) {
    if (all(c("predicted_known_rep", "predicted_missing_rep") %in%
      names(nowcast_summary$nowcast))) {
      nowcast_summary$posterior_nowcast$predicted_all_rep <-
        nowcast_summary$posterior_nowcast$predicted_known_rep +
        nowcast_summary$posterior_nowcast$predicted_missing_rep
    } else if (("predicted_known_rep" %in% names(nowcast_summary$nowcast)) &
      !("predicted_missing_rep" %in% names(nowcast_summary$nowcast))) {
      # assume that there were not missing reference dates, hence "known" and "all" are identical
      nowcast_summary$posterior_nowcast$predicted_all_rep <-
        nowcast_summary$posterior_nowcast$predicted_known_rep
    }
  }

  # Update summary
  nowcast_summary$posterior_nowcast %>%
    group_by(date) %>%
    median_qi(.width = c(0.5, 0.95), .simple_names = F) %>%
    mutate(.width = as.factor(.width)) %>%
    select(-c(.point, .interval)) %>%
    return()
}

#' Get "naive" nowcasts, i.e. just the (downward biased) count of cases observed
#' until now. This is useful for validation, showing how much better the nowcast
#' is compared to a situation in which no nowcasting is done at all (naive
#' nowcast).
get_naive_nowcasts_empirical <- function(ground_truth_df, maxDelay) {
  ground_truth_agg <- ground_truth_df %>%
    group_by(event1_date, event2_date) %>%
    count()

  ground_truth_agg %>%
    group_by(event1_date) %>%
    transmute(nowcast_date = event2_date, nowcast_known_naive = cumsum(n)) %>%
    filter(!is.na(event1_date)) %>%
    ungroup() %>%
    complete(
      event1_date = seq.Date(min(event1_date, na.rm = T),
        max(nowcast_date, na.rm = T),
        by = "1 day"
      ),
      nowcast_date = seq.Date(min(event1_date, na.rm = T),
        max(nowcast_date, na.rm = T),
        by = "1 day"
      )
    ) %>%
    filter(nowcast_date >= event1_date) %>%
    filter(nowcast_date <= event1_date + maxDelay) %>%
    arrange(event1_date, nowcast_date) %>%
    group_by(event1_date) %>%
    fill(nowcast_known_naive, .direction = "downup") %>%
    mutate(nowcast_known_naive = na.fill(nowcast_known_naive, 0)) %>%
    return()
}

#' Add "naive" nowcasts to the simulated results, i.e. just the (downward
#' biased) count of cases observed until now. This is useful for validation,
#' showing how much better the nowcast is compared to a situation in which no
#' nowcasting is done at all (naive nowcast).
#' 
#' @param ground_truth_sim_list A `list` of simulated ground truth values, each
#'   entry in the list represents a different simulation run.
#' @param reference_date The starting reference date to map the simulation time
#'   points to dates.
add_naive_nowcasts_simulated <- function(nowcasts, ground_truth_sim_list, 
                                         reference_date) {
  # note that here we use the maximum delay (not the maximum modeled one), hence
  # the time window, on purpose
  maxDelay <- max(nowcasts$delay, na.rm = T)

  naive_nowcasts <- bind_rows(lapply(ground_truth_sim_list,
                                     function(ground_truth_sim) {
    # known
    ground_truth_agg <- ground_truth_sim$linelist %>%
      filter(onset_known) %>%
      group_by(onset_time, rep_time) %>%
      count()

    known_cases <- ground_truth_agg %>%
      rename(date = onset_time) %>%
      group_by(date) %>%
      transmute(
        date = reference_date + date,
        nowcast_date = reference_date + rep_time,
        nowcast_known_naive = cumsum(n)
      ) %>%
      ungroup() %>%
      complete(
        date = seq.Date(min(date, na.rm = T),
          max(nowcast_date, na.rm = T),
          by = "1 day"
        ),
        nowcast_date = seq.Date(min(date, na.rm = T),
          max(nowcast_date, na.rm = T),
          by = "1 day"
        )
      ) %>%
      filter(nowcast_date >= date) %>%
      filter(nowcast_date <= date + maxDelay) %>%
      arrange(date, nowcast_date) %>%
      group_by(date) %>%
      fill(nowcast_known_naive, .direction = "downup") %>%
      mutate(nowcast_known_naive = na.fill(nowcast_known_naive, 0))

    # all
    ground_truth_agg <- ground_truth_sim$linelist %>%
      group_by(onset_time, rep_time) %>%
      count()

    all_cases <- ground_truth_agg %>%
      rename(date = onset_time) %>%
      group_by(date) %>%
      transmute(
        date = reference_date + date,
        nowcast_date = reference_date + rep_time,
        nowcast_all_naive = cumsum(n)
      ) %>%
      ungroup() %>%
      complete(
        date = seq.Date(min(date, na.rm = T),
          max(nowcast_date, na.rm = T),
          by = "1 day"
        ),
        nowcast_date = seq.Date(min(date, na.rm = T),
          max(nowcast_date, na.rm = T),
          by = "1 day"
        )
      ) %>%
      filter(nowcast_date >= date) %>%
      filter(nowcast_date <= date + maxDelay) %>%
      arrange(date, nowcast_date) %>%
      group_by(date) %>%
      fill(nowcast_all_naive, .direction = "downup") %>%
      mutate(nowcast_all_naive = na.fill(nowcast_all_naive, 0))

    naive_result <- known_cases %>%
      full_join(all_cases, by = c("date", "nowcast_date"))

    return(naive_result)
  }), .id = "dataset_index")

  nowcasts <- nowcasts %>%
    left_join(naive_nowcasts, by = c("nowcast_date", "date", "dataset_index"))
  return(nowcasts)
}

#' Gets a "ground truth" from consolidated, empirical data, and from
#' consolidated nowcasts (i.e. at long longs).
#'
#' It is important to make sure that the data is really consolidated, i.e. all
#' cases have been ' observed for the provided date range
get_ground_truth_empirical <- function(all_nowcasts, empirical_linelist,
                                       maxDelay, consolidation_lags = 14:28) {
  ground_truth_complete <- empirical_linelist %>%
    filter(
      !is.na(event1_date),
      event1_date <= max(all_nowcasts$date, na.rm = T),
      event1_date >= min(all_nowcasts$date, na.rm = T),
      event2_date - event1_date >= 0,
      event2_date - event1_date <= maxDelay
    ) %>%
    count(date = event1_date, name = "nowcast_known_true")

  ground_truth_missing <- empirical_linelist %>%
    filter(
      is.na(event1_date),
      event2_date <= max(all_nowcasts$date, na.rm = T) + maxDelay,
      event2_date >= min(all_nowcasts$date, na.rm = T) - maxDelay
    ) %>%
    count(date = event2_date, name = "missing_rep_true")

  ground_truth_R <- all_nowcasts %>%
    get_consolidated_R_as_true(maxDelay + consolidation_lags)

  ground_truth_cases <- all_nowcasts %>%
    get_consolidated_nowcast_as_true(maxDelay + consolidation_lags)

  ground_truth_all <- list(
    ground_truth_complete,
    ground_truth_missing,
    ground_truth_R,
    ground_truth_cases
  ) %>%
    purrr::reduce(full_join, by = "date")

  return(ground_truth_all)
}

#' Add simulated ground truth values to the nowcast results.
#'
#' @param ground_truth_sim_list A `list` of simulated ground truth values, each
#'   entry in the list represents a different simulation run.
#' @param maxDelay Maximum assumed delay used by the models.
#' @param reference_date The starting reference date to map the simulation time
#'   points to dates.
#' @param full_join Should ground truth values also be added for dates not
#'   present in the nowcasting results?
add_ground_truth_simulated <- function(nowcasts, ground_truth_sim_list, maxDelay,
                                       reference_date, full_join = F) {
  mindate <- min(nowcasts$date, na.rm = T)
  maxdate <- max(nowcasts$date, na.rm = T)

  nowcast_true <- bind_rows(lapply(ground_truth_sim_list, function(ground_truth_sim) {
    data.frame(
      date = reference_date + ground_truth_sim[["process_summary"]]$t,
      nowcast_known_true = ground_truth_sim[["process_summary"]]$onsets_observed,
      nowcast_unknown_true = ground_truth_sim[["process_summary"]]$onsets_hospitalized -
        ground_truth_sim[["process_summary"]]$onsets_observed,
      nowcast_all_true = ground_truth_sim[["process_summary"]]$onsets_hospitalized,
      R_true = ground_truth_sim[["parameters"]]$R
    )
  }), .id = "dataset_index")

  missing_rep_true <- bind_rows(lapply(ground_truth_sim_list, function(ground_truth_sim) {
    ground_truth_sim$linelist %>%
      filter(!onset_known) %>%
      mutate(rep_time = reference_date + rep_time) %>%
      filter(
        rep_time <= maxdate + maxDelay,
        rep_time >= mindate - maxDelay
      ) %>%
      count(rep_time, name = "missing_rep_true")
  }), .id = "dataset_index")

  mean_delay_true <- bind_rows(lapply(ground_truth_sim_list, function(ground_truth_sim) {
    data.frame(
      date = reference_date + ground_truth_sim[["process_summary"]]$t_onset,
      mean_delay_true = dist_get_mean(ground_truth_sim[["parameters"]]$symp_to_rep_dist[ground_truth_sim[["process_summary"]]$t_onset, ])
    )
  }), .id = "dataset_index")

  alpha_true <- bind_rows(lapply(ground_truth_sim_list, function(ground_truth_sim) {
    data.frame(
      date = reference_date + ground_truth_sim[["process_summary"]]$t_onset,
      alpha_true = ground_truth_sim[["parameters"]]$onset_known_prob[ground_truth_sim[["process_summary"]]$t_onset]
    )
  }), .id = "dataset_index")

  if (full_join) {
    join_f <- dplyr::full_join
  } else {
    join_f <- dplyr::left_join
  }

  nowcasts <- nowcasts %>%
    join_f(nowcast_true, by = c("date", "dataset_index")) %>%
    join_f(missing_rep_true, by = c("date" = "rep_time", "dataset_index")) %>%
    join_f(mean_delay_true, by = c("date", "dataset_index")) %>%
    join_f(alpha_true, by = c("date", "dataset_index"))

  return(nowcasts)
}

#' When validating on empirical data, we have no true R. As a proxy, we can use
#' R estimates obtained on the fully consolidated data over all models
#' 
#' @param lags The lags at which to regard data as fully consolidated.
get_consolidated_R_as_true <- function(nowcasts, lags) {
  if ("R" %in% names(nowcasts)) {
    final_estimates <- nowcasts %>%
      filter(as.integer(delay) %in% lags) %>%
      group_by(date) %>%
      summarize(R_true = mean(R, na.rm = T))
    return(final_estimates)
  }
  stop("No R nowcast available")
}

#' Add R estimates obtained on fully consolidated data as proxy ground truth.
#' 
#' @param lags The lags at which to regard data as fully consolidated.
add_consolidated_R_as_true <- function(nowcasts, lags) {
  if ("R" %in% names(nowcasts)) {
    nowcasts <- nowcasts %>%
      left_join(get_consolidated_R_as_true(nowcasts, lags), by = c("date"))
  }
  return(nowcasts)
}

#' When validating on empirical data with missing symptom onset dates, we have
#' no true number of cases. As a proxy, we can use the "nowcast" obtained on the
#' fully consolidated data over all models
#' 
#' @param lags The lags at which to regard data as fully consolidated.
get_consolidated_nowcast_as_true <- function(nowcasts, lags) {
  if ("nowcast_all" %in% names(nowcasts)) {
    final_estimates <- nowcasts %>%
      filter(as.integer(delay) %in% lags) %>%
      group_by(date) %>%
      summarize(nowcast_all_true = mean(nowcast_all, na.rm = T))
    return(final_estimates)
  }
  stop("No nowcast available")
}

#' Add case number estimates obtained on fully consolidated data as proxy ground
#' truth.
#'
#' @param lags The lags at which to regard data as fully consolidated.
add_consolidated_nowcast_as_true <- function(nowcasts, lags) {
  if ("nowcast_all" %in% names(nowcasts)) {
    nowcasts <- nowcasts %>%
      left_join(get_consolidated_nowcast_as_true(nowcasts, lags), by = c("date"))
  }
  return(nowcasts)
}

#' Add various per-observation performance measures to the nowcast results.
#' These are later used as building blocks in `get_metrics()`.
add_performance <- function(nowcasts) {
  vars <- names(nowcasts)

  get_wis_scores <- function(y, dat, alphas = c(0.02, 0.05, seq(0.1, 0.9, 0.1))) {
    wis_matrix <- sapply(alphas, function(alpha) {
      unlist(scoringutils::interval_score(
        y,
        lower = quantile(dat, 0.5 - alpha / 2),
        upper = quantile(dat, 0.5 + alpha / 2),
        interval_range = alpha * 100,
        weigh = T, separate_results = T
      ))
    })
    error <- y - median(dat)
    wis_median <- 0.5 * c(
      "wis" = abs(error),
      "wis_disp" = 0,
      "wis_under" = max(error, 0),
      "wis_over" = max(-error, 0)
    )
    wis_matrix <- cbind(wis_median, wis_matrix)
    return(rowSums(wis_matrix) / (length(alphas) + 0.5))
  }

  # nowcast known
  if ("nowcast_known" %in% vars && "nowcast_known_true" %in% vars) {
    nowcasts$nowcast_known_err <-
      nowcasts$nowcast_known -
      nowcasts$nowcast_known_true
  }
  if ("nowcast_known.upper" %in% vars && "nowcast_known.lower" %in% vars) {
    nowcasts$nowcast_known_range <-
      nowcasts$nowcast_known.upper -
      nowcasts$nowcast_known.lower
  }
  if ("nowcast_known_true" %in% vars && 
      "nowcast_known.upper" %in% vars && 
      "nowcast_known.lower" %in% vars) {
    nowcasts$nowcast_known_within_interval <-
      (nowcasts$nowcast_known_true >= nowcasts$nowcast_known.lower &
        nowcasts$nowcast_known_true <= nowcasts$nowcast_known.upper)
  }
  if ("nowcast_known_true" %in% vars && 
      "nowcast_known_posterior" %in% vars) {
    nowcasts$nowcast_known_crps <- mapply(
      dat = nowcasts$nowcast_known_posterior,
      y = nowcasts$nowcast_known_true, function(dat, y) {
        if (is.na(y) | is.null(dat)) {
          return(NA)
        } else {
          scoringutils::crps_sample(y, matrix(dat, nrow = 1))
        }
      }
    )
    nowcasts$nowcast_known_logS <- mapply(
      dat = nowcasts$nowcast_known_posterior,
      y = nowcasts$nowcast_known_true, function(dat, y) {
        if (is.na(y) | is.null(dat)) {
          return(NA)
        } else {
          scoringutils::logs_sample(y, matrix(dat, nrow = 1))
        }
      }
    )
    nowcast_known_wis <- do.call(rbind, mapply(
      dat = nowcasts$nowcast_known_posterior,
      y = nowcasts$nowcast_known_true, function(dat, y) {
        if (is.na(y) | is.null(dat)) {
          return(NA)
        } else {
          get_wis_scores(y, dat)
        }
      }, SIMPLIFY = F
    ))
    colnames(nowcast_known_wis) <- paste0(
      "nowcast_known_", colnames(nowcast_known_wis)
    )
    nowcasts <- cbind(nowcasts, nowcast_known_wis)
  }

  # nowcast unknown
  if ("nowcast_unknown" %in% vars && "nowcast_unknown_true" %in% vars) {
    nowcasts$nowcast_unknown_err <-
      nowcasts$nowcast_unknown -
      nowcasts$nowcast_unknown_true
  }
  if ("nowcast_unknown.upper" %in% vars && "nowcast_unknown.lower" %in% vars) {
    nowcasts$nowcast_unknown_range <-
      nowcasts$nowcast_unknown.upper -
      nowcasts$nowcast_unknown.lower
  }
  if ("nowcast_unknown_true" %in% vars && 
      "nowcast_unknown.upper" %in% vars && 
      "nowcast_unknown.lower" %in% vars) {
    nowcasts$nowcast_unknown_within_interval <-
      (nowcasts$nowcast_unknown_true >= nowcasts$nowcast_unknown.lower &
        nowcasts$nowcast_unknown_true <= nowcasts$nowcast_unknown.upper)
  }
  if ("nowcast_unknown_true" %in% vars &&
      "nowcast_unknown_posterior" %in% vars) {
    nowcasts$nowcast_unknown_crps <- mapply(
      dat = nowcasts$nowcast_unknown_posterior,
      y = nowcasts$nowcast_unknown_true, function(dat, y) {
        if (is.na(y) | is.null(dat)) {
          return(NA)
        } else {
          crps_sample(y, matrix(dat, nrow = 1))
        }
      }
    )
    nowcasts$nowcast_unknown_logS <- mapply(
      dat = nowcasts$nowcast_unknown_posterior,
      y = nowcasts$nowcast_unknown_true, function(dat, y) {
        if (is.na(y) | is.null(dat)) {
          return(NA)
        } else {
          logs_sample(y, matrix(dat, nrow = 1))
        }
      }
    )
    nowcast_unknown_wis <- do.call(rbind, mapply(
      dat = nowcasts$nowcast_unknown_posterior,
      y = nowcasts$nowcast_unknown_true, function(dat, y) {
        if (is.na(y) | is.null(dat)) {
          return(NA)
        } else {
          get_wis_scores(y, dat)
        }
      }, SIMPLIFY = F
    ))
    colnames(nowcast_unknown_wis) <- paste0(
      "nowcast_unknown_", colnames(nowcast_unknown_wis)
    )
    nowcasts <- cbind(nowcasts, nowcast_unknown_wis)
  }

  # nowcast all
  if ("nowcast_all" %in% vars && "nowcast_all_true" %in% vars) {
    nowcasts$nowcast_all_err <-
      nowcasts$nowcast_all -
      nowcasts$nowcast_all_true
  }
  if ("nowcast_all.upper" %in% vars && "nowcast_all.lower" %in% vars) {
    nowcasts$nowcast_all_range <-
      nowcasts$nowcast_all.upper -
      nowcasts$nowcast_all.lower
  }
  if ("nowcast_all_true" %in% vars &&
      "nowcast_all.upper" %in% vars &&
      "nowcast_all.lower" %in% vars) {
    nowcasts$nowcast_all_within_interval <-
      (nowcasts$nowcast_all_true >= nowcasts$nowcast_all.lower &
        nowcasts$nowcast_all_true <= nowcasts$nowcast_all.upper)
  }
  if ("nowcast_all_true" %in% vars &&
      "nowcast_all_posterior" %in% vars) {
    nowcasts$nowcast_all_crps <- mapply(
      dat = nowcasts$nowcast_all_posterior,
      y = nowcasts$nowcast_all_true, function(dat, y) {
        if (is.na(y) | is.null(dat)) {
          return(NA)
        } else {
          crps_sample(y, matrix(dat, nrow = 1))
        }
      }
    )
    nowcasts$nowcast_all_logS <- mapply(
      dat = nowcasts$nowcast_all_posterior,
      y = nowcasts$nowcast_all_true, function(dat, y) {
        if (is.na(y) | is.null(dat)) {
          return(NA)
        } else {
          logs_sample(y, matrix(dat, nrow = 1))
        }
      }
    )
    nowcast_all_wis <- do.call(rbind, mapply(
      dat = nowcasts$nowcast_all_posterior,
      y = nowcasts$nowcast_all_true, function(dat, y) {
        if (is.na(y) | is.null(dat)) {
          return(NA)
        } else {
          get_wis_scores(y, dat)
        }
      }, SIMPLIFY = F
    ))
    colnames(nowcast_all_wis) <- paste0(
      "nowcast_all_", colnames(nowcast_all_wis)
    )
    nowcasts <- cbind(nowcasts, nowcast_all_wis)
  }

  # missing rep
  if ("predicted_missing_rep" %in% vars &&
      "missing_rep_true" %in% vars) {
    nowcasts$missing_rep_err <-
      nowcasts$predicted_missing_rep -
      nowcasts$missing_rep_true
  }
  if ("predicted_missing_rep.upper" %in% vars &&
      "predicted_missing_rep.lower" %in% vars) {
    nowcasts$missing_rep_range <-
      nowcasts$predicted_missing_rep.upper -
      nowcasts$predicted_missing_rep.lower
  }
  if ("missing_rep_true" %in% vars &&
      "predicted_missing_rep.upper" %in% vars &&
      "predicted_missing_rep.lower" %in% vars) {
    nowcasts$missing_rep_within_interval <-
      (nowcasts$missing_rep_true >= nowcasts$predicted_missing_rep.lower &
        nowcasts$missing_rep_true <= nowcasts$predicted_missing_rep.upper)
  }
  if ("missing_rep_true" %in% vars &&
      "predicted_missing_rep_posterior" %in% vars) {
    nowcasts$missing_rep_crps <- mapply(
      dat = nowcasts$predicted_missing_rep_posterior,
      y = nowcasts$missing_rep_true, function(dat, y) {
        if (is.na(y) | is.null(dat)) {
          return(NA)
        } else {
          crps_sample(y, matrix(dat, nrow = 1))
        }
      }
    )
    nowcasts$missing_rep_logS <- mapply(
      dat = nowcasts$predicted_missing_rep_posterior,
      y = nowcasts$missing_rep_true, function(dat, y) {
        if (is.na(y) | is.null(dat)) {
          return(NA)
        } else {
          logs_sample(y, matrix(dat, nrow = 1))
        }
      }
    )
    missing_rep_wis <- do.call(rbind, mapply(
      dat = nowcasts$predicted_missing_rep_posterior,
      y = nowcasts$missing_rep_true, function(dat, y) {
        if (is.na(y) | is.null(dat)) {
          return(NA)
        } else {
          get_wis_scores(y, dat)
        }
      }, SIMPLIFY = F
    ))
    colnames(missing_rep_wis) <- paste0("missing_rep_", colnames(missing_rep_wis))
    nowcasts <- cbind(nowcasts, missing_rep_wis)
  }

  # mean delay
  if ("mean_delay" %in% vars && "mean_delay_true" %in% vars) {
    nowcasts$mean_delay_err <-
      nowcasts$mean_delay -
      nowcasts$mean_delay_true
  }
  if ("mean_delay.upper" %in% vars && "mean_delay.lower" %in% vars) {
    nowcasts$mean_delay_range <-
      nowcasts$mean_delay.upper -
      nowcasts$mean_delay.lower
  }
  if ("mean_delay_true" %in% vars &&
      "mean_delay.upper" %in% vars &&
      "mean_delay.lower" %in% vars) {
    nowcasts$mean_delay_within_interval <-
      (nowcasts$mean_delay_true >= nowcasts$mean_delay.lower &
        nowcasts$mean_delay_true <= nowcasts$mean_delay.upper)
  }

  # naive nowcast (for known)
  if ("nowcast_known_naive" %in% vars &&
      "nowcast_known_true" %in% vars) {
    nowcasts$nowcast_known_naive_err <-
      nowcasts$nowcast_known_naive -
      nowcasts$nowcast_known_true
  }
  # naive nowcast (for unknown)
  if ("nowcast_unknown_naive" %in% vars &&
      "nowcast_unknown_true" %in% vars) {
    nowcasts$nowcast_unknown_naive_err <-
      nowcasts$nowcast_unknown_naive -
      nowcasts$nowcast_unknown_true
  }
  # naive nowcast (for all)
  if ("nowcast_all_naive" %in% vars &&
      "nowcast_all_true" %in% vars) {
    nowcasts$nowcast_all_naive_err <-
      nowcasts$nowcast_all_naive -
      nowcasts$nowcast_all_true
  }

  # Reproduction number
  if ("R" %in% vars && "R_true" %in% vars) {
    nowcasts$R_err <- nowcasts$R - nowcasts$R_true
  }
  if ("R.upper" %in% vars && "R.lower" %in% vars) {
    nowcasts$R_range <- nowcasts$R.upper - nowcasts$R.lower
  }
  if ("R_true" %in% vars && "R.upper" %in% vars && "R.lower" %in% vars) {
    nowcasts$R_within_interval <-
      (nowcasts$R_true >= nowcasts$R.lower &
        nowcasts$R_true <= nowcasts$R.upper)
  }
  if ("R_true" %in% vars && "R_posterior" %in% vars) {
    nowcasts$R_crps <- mapply(
      dat = nowcasts$R_posterior,
      y = nowcasts$R_true, FUN = function(dat, y) {
        if (is.na(y) | is.null(dat)) {
          return(NA)
        } else {
          crps_sample(y, matrix(dat, nrow = 1))
        }
      }
    )
    nowcasts$R_logS <- mapply(
      dat = nowcasts$R_posterior,
      y = nowcasts$R_true, FUN = function(dat, y) {
        if (is.na(y) | is.null(dat)) {
          return(NA)
        } else {
          logs_sample(y, matrix(dat, nrow = 1))
        }
      }
    )
    R_wis <- do.call(rbind, mapply(
      dat = nowcasts$R_posterior,
      y = nowcasts$R_true, FUN = function(dat, y) {
        if (is.na(y) | is.null(dat)) {
          return(NA)
        } else {
          get_wis_scores(y, dat)
        }
      }, SIMPLIFY = FALSE
    ))
    colnames(R_wis) <- paste0("R_", colnames(R_wis))
    nowcasts <- cbind(nowcasts, R_wis)
  }

  return(nowcasts)
}

#' Helper function to temporarily add missing performance columns as NA
add_missing_performance <- function(nowcasts) {
  columns <- c(
    "nowcast_known",
    "nowcast_unknown",
    "nowcast_all",
    "missing_rep",
    "mean_delay",
    "R"
  )
  suffixes <- c(
    "",
    "_true",
    "_err",
    "_naive_err",
    "_range",
    "_within_interval",
    "_point_viol",
    "_any_viol",
    "_crps",
    "_logS",
    "_wis",
    "_wis_disp",
    "_wis_under",
    "_wis_over",
    "_posterior"
  )
  perfcolnames <- paste0(
    rep(columns, each = length(suffixes)),
    rep(suffixes, length(columns))
  )
  cols <- setNames(rep(NA, length(perfcolnames)), perfcolnames)

  nowcasts <- add_column(
    nowcasts, !!!cols[setdiff(names(cols), names(nowcasts))]
  )
  return(nowcasts)
}

#' Compute various metrics for evaluation of point and interval forecasts
get_metrics <- function(nowcasts) {
  nowcasts <- add_missing_performance(nowcasts)

  metrics <- nowcasts %>%
    summarize(
      # nowcast known
      nowcast_known_MAE = mean(abs(nowcast_known_err), na.rm = T),
      nowcast_known_naive_MAE = mean(abs(nowcast_known_naive_err), na.rm = T),
      nowcast_known_MASE = nowcast_known_MAE / nowcast_known_naive_MAE,
      nowcast_known_MAPE = mean(abs(nowcast_known_err) / abs(nowcast_known_true), na.rm = T),
      nowcast_known_mean_range = mean(nowcast_known_range, na.rm = T),
      nowcast_known_coverage = sum(nowcast_known_within_interval, na.rm = T) / sum(!is.na(nowcast_known_within_interval)),
      nowcast_known_point_consistency = mean(1 - nowcast_known_point_viol, na.rm = T),
      nowcast_known_point_consistency_full = 1 - sum(nowcast_known_point_viol > 0, na.rm = T) / sum(!is.na(nowcast_known_point_viol)),
      nowcast_known_interval_consistency = mean(1 - nowcast_known_any_viol, na.rm = T),
      nowcast_known_interval_consistency_full = 1 - sum(nowcast_known_any_viol > 0, na.rm = T) / sum(!is.na(nowcast_known_any_viol)),
      nowcast_known_mean_crps = mean(nowcast_known_crps, na.rm = T),
      nowcast_known_mean_crps_scaled = mean(abs(nowcast_known_crps) / abs(nowcast_known_true), na.rm = T),
      nowcast_known_mean_logS = mean(nowcast_known_logS, na.rm = T),
      nowcast_known_mean_wis = mean(nowcast_known_wis, na.rm = T),
      nowcast_known_mean_wis_disp = mean(nowcast_known_wis_disp, na.rm = T),
      nowcast_known_mean_wis_under = mean(nowcast_known_wis_under, na.rm = T),
      nowcast_known_mean_wis_over = mean(nowcast_known_wis_over, na.rm = T),
      # nowcast unknown
      nowcast_unknown_MAE = mean(abs(nowcast_unknown_err), na.rm = T),
      nowcast_unknown_naive_MAE = mean(abs(nowcast_unknown_naive_err), na.rm = T),
      nowcast_unknown_MASE = nowcast_unknown_MAE / nowcast_unknown_naive_MAE,
      nowcast_unknown_MAPE = mean(abs(nowcast_unknown_err) / abs(nowcast_unknown_true), na.rm = T),
      nowcast_unknown_mean_range = mean(nowcast_unknown_range, na.rm = T),
      nowcast_unknown_coverage = sum(nowcast_unknown_within_interval, na.rm = T) / sum(!is.na(nowcast_unknown_within_interval)),
      nowcast_unknown_point_consistency = mean(1 - nowcast_unknown_point_viol, na.rm = T),
      nowcast_unknown_point_consistency_full = 1 - sum(nowcast_unknown_point_viol > 0, na.rm = T) / sum(!is.na(nowcast_unknown_point_viol)),
      nowcast_unknown_interval_consistency = mean(1 - nowcast_unknown_any_viol, na.rm = T),
      nowcast_unknown_interval_consistency_full = 1 - sum(nowcast_unknown_any_viol > 0, na.rm = T) / sum(!is.na(nowcast_unknown_any_viol)),
      nowcast_unknown_mean_crps = mean(nowcast_unknown_crps, na.rm = T),
      nowcast_unknown_mean_crps_scaled = mean(abs(nowcast_unknown_crps) / abs(nowcast_unknown_true), na.rm = T),
      nowcast_unknown_mean_logS = mean(nowcast_unknown_logS, na.rm = T),
      nowcast_unknown_mean_wis = mean(nowcast_unknown_wis, na.rm = T),
      nowcast_unknown_mean_wis_disp = mean(nowcast_unknown_wis_disp, na.rm = T),
      nowcast_unknown_mean_wis_under = mean(nowcast_unknown_wis_under, na.rm = T),
      nowcast_unknown_mean_wis_over = mean(nowcast_unknown_wis_over, na.rm = T),
      # nowcast all
      nowcast_all_MAE = mean(abs(nowcast_all_err), na.rm = T),
      nowcast_all_naive_MAE = mean(abs(nowcast_all_naive_err), na.rm = T),
      nowcast_all_MASE = nowcast_all_MAE / nowcast_all_naive_MAE,
      nowcast_all_MAPE = mean(abs(nowcast_all_err) / abs(nowcast_all_true), na.rm = T),
      nowcast_all_mean_range = mean(nowcast_all_range, na.rm = T),
      nowcast_all_coverage = sum(nowcast_all_within_interval, na.rm = T) / sum(!is.na(nowcast_all_within_interval)),
      nowcast_all_point_consistency = mean(1 - nowcast_all_point_viol, na.rm = T),
      nowcast_all_point_consistency_full = 1 - sum(nowcast_all_point_viol > 0, na.rm = T) / sum(!is.na(nowcast_all_point_viol)),
      nowcast_all_interval_consistency = mean(1 - nowcast_all_any_viol, na.rm = T),
      nowcast_all_interval_consistency_full = 1 - sum(nowcast_all_any_viol > 0, na.rm = T) / sum(!is.na(nowcast_all_any_viol)),
      nowcast_all_mean_crps = mean(nowcast_all_crps, na.rm = T),
      nowcast_all_mean_crps_scaled = mean(abs(nowcast_all_crps) / abs(nowcast_all_true), na.rm = T),
      nowcast_all_mean_logS = mean(nowcast_all_logS, na.rm = T),
      nowcast_all_mean_wis = mean(nowcast_all_wis, na.rm = T),
      nowcast_all_mean_wis_disp = mean(nowcast_all_wis_disp, na.rm = T),
      nowcast_all_mean_wis_under = mean(nowcast_all_wis_under, na.rm = T),
      nowcast_all_mean_wis_over = mean(nowcast_all_wis_over, na.rm = T),
      # missing rep
      missing_rep_MAE = mean(abs(missing_rep_err), na.rm = T),
      missing_rep_MAPE = mean(abs(missing_rep_err) / abs(missing_rep_true), na.rm = T),
      missing_rep_mean_range = mean(missing_rep_range, na.rm = T),
      missing_rep_coverage = sum(missing_rep_within_interval, na.rm = T) / sum(!is.na(missing_rep_within_interval)),
      missing_rep_mean_crps = mean(missing_rep_crps, na.rm = T),
      missing_rep_mean_crps_scaled = mean(abs(missing_rep_crps) / abs(missing_rep_true), na.rm = T),
      missing_rep_mean_logS = mean(missing_rep_logS, na.rm = T),
      missing_rep_mean_wis = mean(missing_rep_wis, na.rm = T),
      missing_rep_mean_wis_disp = mean(missing_rep_wis_disp, na.rm = T),
      missing_rep_mean_wis_under = mean(missing_rep_wis_under, na.rm = T),
      missing_rep_mean_wis_over = mean(missing_rep_wis_over, na.rm = T),
      # mean delay
      mean_delay_MAE = mean(abs(mean_delay_err), na.rm = T),
      mean_delay_MAPE = mean(abs(mean_delay_err) / abs(mean_delay_true), na.rm = T),
      mean_delay_mean_range = mean(mean_delay_range, na.rm = T),
      mean_delay_coverage = sum(mean_delay_within_interval, na.rm = T) / sum(!is.na(mean_delay_within_interval)),
      # R
      R_MAE = mean(abs(R_err), na.rm = T),
      R_MAPE = mean(abs(R_err) / abs(R_true), na.rm = T),
      R_mean_range = mean(R_range, na.rm = T),
      R_coverage = sum(R_within_interval, na.rm = T) / sum(!is.na(R_within_interval)),
      R_above_one_correct = sum((R_true > 1 & R > 1) | (R_true < 1 & R < 1), na.rm = T) / sum(!is.na(R_true) & !is.na(R)),
      R_point_consistency = mean(1 - R_point_viol, na.rm = T),
      R_point_consistency_full = 1 - sum(R_point_viol > 0, na.rm = T) / sum(!is.na(R_point_viol)),
      R_interval_consistency = mean(1 - R_any_viol, na.rm = T),
      R_interval_consistency_full = 1 - sum(R_any_viol > 0, na.rm = T) / sum(!is.na(R_any_viol)),
      R_mean_crps = mean(R_crps, na.rm = T),
      R_mean_crps_scaled = mean(abs(R_crps) / abs(R_true), na.rm = T),
      R_mean_logS = mean(R_logS, na.rm = T),
      R_mean_wis = mean(R_wis, na.rm = T),
      R_mean_wis_disp = mean(R_wis_disp, na.rm = T),
      R_mean_wis_under = mean(R_wis_under, na.rm = T),
      R_mean_wis_over = mean(R_wis_over, na.rm = T),
      .groups = "drop"
    )
  metrics <- metrics %>% select(where(~ !all(is.na(.x))))

  return(metrics)
}

#' Compute various rolling metrics for evaluation of point and interval forecasts
#'
#' @param window_size Number of consecutive days to include in each window
get_metrics_rolling <- function(nowcasts, window_size = 7) {
  nowcasts <- add_missing_performance(nowcasts)

  metrics <- nowcasts %>%
    filter(!is.na(delay)) %>%
    arrange(date) %>%
    summarize(
      date_from = rollmaxr(date, window_size) - window_size,
      date_to = rollmaxr(date, window_size),
      # nowcast known
      nowcast_known_MAE = rollmeanr(abs(nowcast_known_err), window_size, na.rm = T),
      nowcast_known_naive_MAE = rollmeanr(abs(nowcast_known_naive_err), window_size, na.rm = T),
      nowcast_known_MASE = nowcast_known_MAE / nowcast_known_naive_MAE,
      nowcast_known_MAPE = rollmeanr(abs(nowcast_known_err) / abs(nowcast_known_true), window_size, na.rm = T),
      nowcast_known_mean_range = rollmeanr(nowcast_known_range, window_size, na.rm = T),
      nowcast_known_coverage = rollsumr(nowcast_known_within_interval, window_size, na.rm = T) / rollsumr(!is.na(nowcast_known_within_interval), window_size),
      nowcast_known_point_consistency = rollmeanr(1 - nowcast_known_point_viol, window_size, na.rm = T),
      nowcast_known_point_consistency_full = 1 - rollsumr(nowcast_known_point_viol > 0, window_size, na.rm = T) / rollsumr(!is.na(nowcast_known_point_viol), window_size),
      nowcast_known_interval_consistency = rollmeanr(1 - nowcast_known_any_viol, window_size, na.rm = T),
      nowcast_known_interval_consistency_full = 1 - rollsumr(nowcast_known_any_viol > 0, window_size, na.rm = T) / rollsumr(!is.na(nowcast_known_any_viol), window_size),
      nowcast_known_mean_crps = rollmeanr(nowcast_known_crps, window_size, na.rm = T),
      nowcast_known_mean_logS = rollmeanr(nowcast_known_logS, window_size, na.rm = T),
      # nowcast unknown
      nowcast_unknown_MAE = rollmeanr(abs(nowcast_unknown_err), window_size, na.rm = T),
      nowcast_unknown_naive_MAE = rollmeanr(abs(nowcast_unknown_naive_err), window_size, na.rm = T),
      nowcast_unknown_MASE = nowcast_unknown_MAE / nowcast_unknown_naive_MAE,
      nowcast_unknown_MAPE = rollmeanr(abs(nowcast_unknown_err) / abs(nowcast_unknown_true), window_size, na.rm = T),
      nowcast_unknown_mean_range = rollmeanr(nowcast_unknown_range, window_size, na.rm = T),
      nowcast_unknown_coverage = rollsumr(nowcast_unknown_within_interval, window_size, na.rm = T) / rollsumr(!is.na(nowcast_unknown_within_interval), window_size),
      nowcast_unknown_point_consistency = rollmeanr(1 - nowcast_unknown_point_viol, window_size, na.rm = T),
      nowcast_unknown_point_consistency_full = 1 - rollsumr(nowcast_unknown_point_viol > 0, window_size, na.rm = T) / rollsumr(!is.na(nowcast_unknown_point_viol), window_size),
      nowcast_unknown_interval_consistency = rollmeanr(1 - nowcast_unknown_any_viol, window_size, na.rm = T),
      nowcast_unknown_interval_consistency_full = 1 - rollsumr(nowcast_unknown_any_viol > 0, window_size, na.rm = T) / rollsumr(!is.na(nowcast_unknown_any_viol), window_size),
      nowcast_unknown_mean_crps = rollmeanr(nowcast_unknown_crps, window_size, na.rm = T),
      nowcast_unknown_mean_logS = rollmeanr(nowcast_unknown_logS, window_size, na.rm = T),
      # nowcast all
      nowcast_all_MAE = rollmeanr(abs(nowcast_all_err), window_size, na.rm = T),
      nowcast_all_naive_MAE = rollmeanr(abs(nowcast_all_naive_err), window_size, na.rm = T),
      nowcast_all_MASE = nowcast_all_MAE / nowcast_all_naive_MAE,
      nowcast_all_MAPE = rollmeanr(abs(nowcast_all_err) / abs(nowcast_all_true), window_size, na.rm = T),
      nowcast_all_mean_range = rollmeanr(nowcast_all_range, window_size, na.rm = T),
      nowcast_all_coverage = rollsumr(nowcast_all_within_interval, window_size, na.rm = T) / rollsumr(!is.na(nowcast_all_within_interval), window_size),
      nowcast_all_point_consistency = rollmeanr(1 - nowcast_all_point_viol, window_size, na.rm = T),
      nowcast_all_point_consistency_full = 1 - rollsumr(nowcast_all_point_viol > 0, window_size, na.rm = T) / rollsumr(!is.na(nowcast_all_point_viol), window_size),
      nowcast_all_interval_consistency = rollmeanr(1 - nowcast_all_any_viol, window_size, na.rm = T),
      nowcast_all_interval_consistency_full = 1 - rollsumr(nowcast_all_any_viol > 0, window_size, na.rm = T) / rollsumr(!is.na(nowcast_all_any_viol), window_size),
      nowcast_all_mean_crps = rollmeanr(nowcast_all_crps, window_size, na.rm = T),
      nowcast_all_mean_logS = rollmeanr(nowcast_all_logS, window_size, na.rm = T),
      # missing rep
      missing_rep_MAE = rollmeanr(abs(missing_rep_err), window_size, na.rm = T),
      missing_rep_MAPE = rollmeanr(abs(missing_rep_err) / abs(missing_rep_true), window_size, na.rm = T),
      missing_rep_mean_range = rollmeanr(missing_rep_range, window_size, na.rm = T),
      missing_rep_coverage = rollsumr(missing_rep_within_interval, window_size, na.rm = T) / rollsumr(!is.na(missing_rep_within_interval), window_size),
      missing_rep_mean_crps = rollmeanr(missing_rep_crps, window_size, na.rm = T),
      missing_rep_mean_logS = rollmeanr(missing_rep_logS, window_size, na.rm = T),
      # mean delay
      mean_delay_MAE = rollmeanr(abs(mean_delay_err), window_size, na.rm = T),
      mean_delay_MAPE = rollmeanr(abs(mean_delay_err) / abs(mean_delay_true), window_size, na.rm = T),
      mean_delay_mean_range = rollmeanr(mean_delay_range, window_size, na.rm = T),
      mean_delay_coverage = rollsumr(mean_delay_within_interval, window_size, na.rm = T) / rollsumr(!is.na(mean_delay_within_interval), window_size),
      # R
      R_MAE = rollmeanr(abs(R_err), window_size, na.rm = T),
      R_MAPE = rollmeanr(abs(R_err) / abs(R_true), window_size, na.rm = T),
      R_mean_range = rollmeanr(R_range, window_size, na.rm = T),
      R_coverage = rollsumr(R_within_interval, window_size, na.rm = T) / rollsumr(!is.na(R_within_interval), window_size),
      R_point_consistency = rollmeanr(1 - R_point_viol, window_size, na.rm = T),
      R_point_consistency_full = 1 - rollsumr(R_point_viol > 0, window_size, na.rm = T) / rollsumr(!is.na(R_point_viol), window_size),
      R_interval_consistency = rollmeanr(1 - R_any_viol, window_size, na.rm = T),
      R_interval_consistency_full = 1 - rollsumr(R_any_viol > 0, window_size, na.rm = T) / rollsumr(!is.na(R_any_viol), window_size),
      R_mean_crps = rollmeanr(R_crps, window_size, na.rm = T),
      R_mean_logS = rollmeanr(R_logS, window_size, na.rm = T),
      .groups = "drop"
    )

  metrics <- metrics %>% select(where(~ !all(is.na(.x))))

  return(metrics)
}

#' Add measures of nowcast consistency over time
#' For example: is a point nowcast for a given date within the interval of earlier nowcasts for that date?
add_consistency <- function(nowcasts, slack = 2, R_slack = 0.1) {
  if ("nowcast_known" %in% names(nowcasts)) {
    nowcast_known_interval_viols <- nowcasts %>%
      select(dataset_index, date, delay, nowcast_known, nowcast_known.lower, nowcast_known.upper) %>%
      inner_join(nowcasts %>% select(dataset_index, date, delay, nowcast_known, nowcast_known.lower, nowcast_known.upper),
        by = c("dataset_index" = "dataset_index", "date" = "date"), suffix = c("", ".temp")
      ) %>%
      filter(delay < delay.temp, delay < maxDelay) %>%
      group_by(dataset_index, date, delay) %>%
      summarize(
        nowcast_known_upper_viol = sum(nowcast_known.upper.temp > (nowcast_known.upper + slack), na.rm = T) / sum(!is.na(nowcast_known.upper.temp) & !is.na(nowcast_known.upper)), # How often is the upper interval bound of later nowcasts above the current upper bound
        nowcast_known_lower_viol = sum(nowcast_known.lower.temp < (nowcast_known.lower - slack), na.rm = T) / sum(!is.na(nowcast_known.lower.temp) & !is.na(nowcast_known.lower)), # How often is the lower interval bound of later nowcasts below the current lower bound
        nowcast_known_any_viol = sum((nowcast_known.upper.temp > nowcast_known.upper + slack) | (nowcast_known.lower.temp < nowcast_known.lower - slack), na.rm = T) / sum(!is.na(nowcast_known.upper.temp) & !is.na(nowcast_known.upper) & !is.na(nowcast_known.lower.temp) & !is.na(nowcast_known.lower)), # How often is either the upper or lower bound violated
        nowcast_known_point_viol = sum(nowcast_known.temp > (nowcast_known.upper + slack) | nowcast_known.temp < (nowcast_known.lower - slack), na.rm = T) / sum(!is.na(nowcast_known.temp) & !is.na(nowcast_known.upper) & !is.na(nowcast_known.lower)), # How often are the point predictions of later nowcasts outside of the current uncertainty interval
        .groups = "drop"
      )

    nowcasts <- nowcasts %>% left_join(nowcast_known_interval_viols, by = c("dataset_index", "date", "delay"))
  }

  if ("nowcast_unknown" %in% names(nowcasts)) {
    nowcast_unknown_interval_viols <- nowcasts %>%
      select(dataset_index, date, delay, nowcast_unknown, nowcast_unknown.lower, nowcast_unknown.upper) %>%
      inner_join(nowcasts %>% select(dataset_index, date, delay, nowcast_unknown, nowcast_unknown.lower, nowcast_unknown.upper),
        by = c("dataset_index" = "dataset_index", "date" = "date"), suffix = c("", ".temp")
      ) %>%
      filter(delay < delay.temp, delay < maxDelay) %>%
      group_by(dataset_index, date, delay) %>%
      summarize(
        nowcast_unknown_upper_viol = sum(nowcast_unknown.upper.temp > (nowcast_unknown.upper + slack), na.rm = T) / sum(!is.na(nowcast_unknown.upper.temp) & !is.na(nowcast_unknown.upper)),
        nowcast_unknown_lower_viol = sum(nowcast_unknown.lower.temp < (nowcast_unknown.lower - slack), na.rm = T) / sum(!is.na(nowcast_unknown.lower.temp) & !is.na(nowcast_unknown.lower)),
        nowcast_unknown_any_viol = sum((nowcast_unknown.upper.temp > nowcast_unknown.upper + slack) | (nowcast_unknown.lower.temp < nowcast_unknown.lower - slack), na.rm = T) / sum(!is.na(nowcast_unknown.upper.temp) & !is.na(nowcast_unknown.upper) & !is.na(nowcast_unknown.lower.temp) & !is.na(nowcast_unknown.lower)),
        nowcast_unknown_point_viol = sum(nowcast_unknown.temp > (nowcast_unknown.upper + slack) | nowcast_unknown.temp < (nowcast_unknown.lower - slack), na.rm = T) / sum(!is.na(nowcast_unknown.temp) & !is.na(nowcast_unknown.upper) & !is.na(nowcast_unknown.lower)),
        .groups = "drop"
      )

    nowcasts <- nowcasts %>% left_join(nowcast_unknown_interval_viols, by = c("dataset_index", "date", "delay"))
  }

  if ("R" %in% names(nowcasts)) {
    R_interval_viols <- nowcasts %>%
      select(dataset_index, date, delay, R, R.lower, R.upper) %>%
      inner_join(nowcasts %>% select(dataset_index, date, delay, R, R.lower, R.upper),
        by = c("dataset_index" = "dataset_index", "date" = "date"), suffix = c("", ".temp")
      ) %>%
      filter(delay < delay.temp, delay < maxDelay) %>%
      group_by(dataset_index, date, delay) %>%
      summarize(
        R_upper_viol = sum(R.upper.temp > (R.upper + R_slack), na.rm = T) / sum(!is.na(R.upper.temp) & !is.na(R.upper)),
        R_lower_viol = sum(R.lower.temp < (R.lower - R_slack), na.rm = T) / sum(!is.na(R.lower.temp) & !is.na(R.lower)),
        R_any_viol = sum((R.upper.temp > R.upper + R_slack) | (R.lower.temp < R.lower - R_slack), na.rm = T) / sum(!is.na(R.upper.temp) & !is.na(R.upper) & !is.na(R.lower.temp) & !is.na(R.lower)),
        R_point_viol = sum(R.temp > (R.upper + R_slack) | R.temp < (R.lower - R_slack), na.rm = T) / sum(!is.na(R.temp) & !is.na(R.upper) & !is.na(R.lower)),
        .groups = "drop"
      )

    nowcasts <- nowcasts %>% left_join(R_interval_viols, by = c("dataset_index", "date", "delay"))
  }

  return(nowcasts)
}

# Plotting

# General performance plot across models, phases, for selected lags
plot_performance_select <- function(m_res_list, outcome_name, models, metric, metric_name, delays_select) {
  bind_rows(m_res_list, .id = "model") %>%
    mutate(model = paste(model, R_model)) %>%
    filter(model %in% models, nowcast_timing_list == paste0(delays_select, collapse = ",")) %>%
    mutate(model = recode_factor(model, !!!model_names, .ordered = T)) %>%
    select(dataset_index, model, phase, starts_with(outcome_name)) %>%
    select(-contains("naive")) %>%
    select(dataset_index, model, phase, contains(c("crps", "logS", "coverage", "MAE", "wis"))) %>%
    pivot_longer(-c(dataset_index, model, phase)) %>%
    mutate(outcome = str_extract(name, outcome_name), name = str_remove(name, paste0(outcome_name, "_"))) %>%
    filter(name == metric, !is.na(value)) %>%
    ggplot(aes(y = value, color = model)) +
    geom_boxplot(outlier.shape = NA) +
    facet_wrap(~phase, nrow = 1, scales = "fixed") +
    theme_bw() +
    theme(
      legend.position = "none",
      strip.background = element_rect(fill = NA),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
    ) +
    ylab(metric_name) +
    coord_cartesian(ylim = c(0, 0.5)) +
    scale_x_continuous(breaks = NULL) +
    scale_color_manual(
      name = "Approach", values = model_colors[names(model_colors) %in% model_names[models]]
    ) %>%
    return()
}


#' WIS performance plot across models, phases, for selected lags
plot_performance_wis_select <- function(m_res_list, outcome_name, models, delays_select) {
  plotdat <- bind_rows(m_res_list, .id = "model") %>%
    mutate(model = paste(model, R_model)) %>%
    filter(
      model %in% models,
      nowcast_timing_list == paste0(delays_select, collapse = ",")
    ) %>%
    mutate(model = recode_factor(model, !!!model_names, .ordered = T)) %>%
    select(dataset_index, model, phase, starts_with(outcome_name)) %>%
    select(-contains("naive")) %>%
    select(dataset_index, model, phase, contains(c("wis"))) %>%
    pivot_longer(-c(dataset_index, model, phase)) %>%
    mutate(
      outcome = str_extract(name, outcome_name),
      name = str_remove(name, paste0(outcome_name, "_"))
    ) %>%
    mutate(name = factor(name, ordered = T, levels = c(
      "mean_wis", "mean_wis_over", "mean_wis_disp", "mean_wis_under"
    )))
  
  plotdat %>%
    filter(name != "mean_wis") %>%
    group_by(model, phase, name) %>%
    summarize(mean_value = mean(value), .groups = "keep") %>%
    filter(!is.na(mean_value)) %>%
    mutate(model_x = as.integer(model)) %>%
    ggplot(aes(x = model_x)) +
    geom_col_pattern(aes(y = mean_value, pattern_color = model, pattern = name),
                     fill = "white", colour = "black", position = "stack"
    ) +
    facet_wrap(~phase, nrow = 1, scales = "free_y") +
    theme_bw() +
    theme(
      legend.position = "none",
      strip.background = element_rect(fill = NA),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank()
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0, ifelse(FALSE, 0.5, 0.05)))
    ) +
    ylab(ifelse(outcome_name == "nowcast_all",
                expression(bar(WIS)(hat(N)[t])),
                expression(bar(WIS)(hat(R)[t]))
    )) +
    scale_pattern_discrete(choices = c("stripe", "circle", "crosshatch")) +
    scale_pattern_color_manual(
      name = "Approach",
      values = model_colors[names(model_colors) %in% model_names[models]]
    ) %>%
    return()
}

plot_performance_wis_best <- function(m_res_list, outcome_name, models, delays_select) {
  plotdat <- bind_rows(m_res_list, .id = "model") %>%
    mutate(model = paste(model, R_model)) %>%
    filter(
      model %in% models,
      nowcast_timing_list == paste0(delays_select, collapse = ",")
    ) %>%
    mutate(model = recode_factor(model, !!!model_names, .ordered = T)) %>%
    select(dataset_index, model, phase, starts_with(outcome_name)) %>%
    select(-contains("naive")) %>%
    select(dataset_index, model, phase, contains(c("wis"))) %>%
    pivot_longer(-c(dataset_index, model, phase)) %>%
    mutate(
      outcome = str_extract(name, outcome_name),
      name = str_remove(name, paste0(outcome_name, "_"))
    ) %>%
    mutate(name = factor(name, ordered = T, levels = c(
      "mean_wis", "mean_wis_over", "mean_wis_disp", "mean_wis_under"
    )))
  
  best_percentage <- plotdat %>%
    filter(name == "mean_wis") %>%
    ungroup() %>%
    group_by(dataset_index, phase) %>%
    mutate(best = value == min(value, na.rm = T)) %>%
    group_by(model, phase, name) %>%
    summarize(
      mean_value = mean(value, na.rm = T),
      percent_best = sum(best, na.rm = T) / sum(!is.na(best)),
      .groups = "keep"
    ) %>%
    group_by(phase, name) %>%
    # normalize to account for double-counting
    mutate(percent_best = round(100 * percent_best / sum(percent_best))) %>%
    mutate(model_x = as.integer(model))
  
  best_percentage %>%
    ggplot() +
    geom_bar(aes(x = percent_best, y = 0, fill = forcats::fct_rev(model)),
             alpha = 0.8,
             orientation = "y",
             position = "fill",
             stat = "identity"
    ) +
    facet_wrap(~phase, nrow = 1, scales = "free_y") +
    annotation_custom(
      grid::textGrob(
        label = "Percent best", gp = grid::gpar(fontsize = 9), vjust = 1.5
      ),
      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -Inf
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      strip.background = element_rect(fill = NA),
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank()
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0))) +
    scale_x_continuous(expand = expansion(mult = c(0, 0))) +
    coord_cartesian(clip = "off") +
    scale_fill_manual(
      name = "Approach",
      values = model_colors[names(model_colors) %in% model_names[models]]
    ) %>%
    return()
}

plot_validation_select <- function(dataset_id,
                                   models,
                                   R_models = c("renewal", "epiestim"),
                                   feather_days,
                                   feather_maxdelay,
                                   feather_margin,
                                   lags_list = list(
                                     "Nowcasts in same week" = 0:6,
                                     "Nowcasts one week after" = 7:13,
                                     "Nowcasts two weeks after" = 14:20
                                   ),
                                   plot_type = "Cases",
                                   performance_metric = "wis",
                                   date_breaks = waiver(),
                                   has_missing_onsets = NULL,
                                   percent_best = T,
                                   truth_data_index = 1,
                                   categories_name = "Approach",
                                   date_labels = "%b %d",
                                   ground_truth_linetype = "solid",
                                   legends.rel_widths = NULL,
                                   legend_plot.rel_heights = c(1, 9.5)) {
  intermediate_plots <- list()
  
  minimum_x <- reference_date + feather_days[1] - feather_maxdelay - 1
  maximum_x <- reference_date + feather_days[1] + feather_margin
  
  if (is.null(legends.rel_widths)) {
    if (performance_metric == "wis") {
      legends.rel_widths <- c(0.5, 0.5)
    } else if (has_missing_onsets & (plot_type != "R")) {
      legends.rel_widths <- c(0.35, 0.45, 0.2)
    } else {
      legends.rel_widths <- c(0.25, 0.5, 0.25)
    }
  }
  
  
  for (i in 1:length(lags_list)) {
    # data preparation
    
    lags <- lags_list[[i]]
    current_nowcast_dates <- phases %>%
      filter(nowcast_timing_list == paste0(lags, collapse = ",")) %>%
      distinct(nowcast_date) %>%
      pull()
    
    
    nowcast_data <- bind_rows(lapply(n_res_list[models], function(df) {
      df %>%
        filter(
          dataset_index == dataset_id,
          R_model %in% R_models,
          nowcast_date %in% (current_nowcast_dates)
        )
    }), .id = "model") %>%
      mutate(
        model = paste(model, R_model),
        group = paste(model, nowcast_date),
        .after = model
      )
    
    models_select <- unique(nowcast_data$model)
    if (is.null(has_missing_onsets)) {
      has_missing_onsets <- any(str_detect(models_select, "miss"))
    }
    
    if ("nowcast_all_naive" %in% names(nowcast_data)) {
      naive_nowcast_data <- nowcast_data %>%
        distinct(nowcast_date, date, nowcast_all_naive)
    }
    
    if (plot_type == "R") {
      onset_legend <- R_legend
    } else {
      if (has_missing_onsets) {
        naive_nowcast_known_data <- nowcast_data %>%
          distinct(nowcast_date, date, nowcast_known_naive)
        onset_legend <- symptom_onset_legend_with_missing
      } else if ("nowcast_all_naive" %in% names(nowcast_data)) {
        naive_nowcast_known_data <- naive_nowcast_data %>%
          rename(nowcast_known_naive = nowcast_all_naive)
        onset_legend <- symptom_onset_legend
      } else {
        onset_legend <- symptom_onset_legend # might want to change this
      }
    }
    
    truth_data <- n_res_list[[ifelse(is.null(truth_data_index), 1, truth_data_index)]] %>%
      filter(dataset_index == dataset_id) %>%
      group_by(date) %>%
      summarize(across(ends_with("true"), first)) %>%
      uncount(length(unique(nowcast_data$nowcast_date))) %>%
      mutate(nowcast_date = rep(unique(nowcast_data$nowcast_date), length.out = n()))
    
    add_feather_date <- function(df) {
      nowcast_date_renaming_feather <- setNames(
        as.character(reference_date + feather_days),
        as.character(current_nowcast_dates)
      )
      df <- df %>%
        mutate(feather_date = as.Date(recode(
          as.character(nowcast_date),
          !!!nowcast_date_renaming_feather
        )))
      return(df)
    }
    
    recode_for_plot <- function(df) {
      nowcast_date_renaming_lag <- setNames(
        names(feather_days),
        as.character(current_nowcast_dates)
      )
      if ("model" %in% names(df)) {
        df <- df %>%
          mutate(model = recode_factor(model, !!!model_names, .ordered = T))
      }
      df <- df %>% mutate(
        nowcast_date = recode_factor(as.character(nowcast_date),
                                     !!!nowcast_date_renaming_lag,
                                     .ordered = T
        ),
      )
      return(df)
    }
    
    nowcast_data_cases <- nowcast_data %>%
      add_feather_date() %>%
      filter(delay <= feather_maxdelay + lags[1] + 1) %>%
      recode_for_plot()
    if (exists("truth_data")) {
      truth_data_cases <- truth_data %>%
        add_feather_date() %>%
        filter(
          date <= feather_date + feather_margin,
          date >= feather_date - feather_maxdelay
        ) %>%
        recode_for_plot()
    }
    if (exists("naive_nowcast_data")) {
      naive_nowcast_data <- naive_nowcast_data %>%
        add_feather_date() %>%
        filter(date >= feather_date - feather_maxdelay) %>%
        recode_for_plot()
    }
    if (exists("naive_nowcast_known_data")) {
      naive_nowcast_known_data <- naive_nowcast_known_data %>%
        add_feather_date() %>%
        filter(date >= feather_date - feather_maxdelay) %>%
        recode_for_plot()
    }
    
    # plot nowcasts
    if (plot_type == "Cases") {
      minimum_y <- -1
      maximum_y <- min(
        max(nowcast_data_cases$nowcast_all.upper, na.rm = T),
        max(c(0, truth_data_cases$nowcast_all_true, na.rm = T)) * 1.3,
        na.rm = T
      ) + 10
      
      PlotNowcast <- nowcast_data_cases %>%
        {
          ggplot(data = .) +
            {
              if (exists("naive_nowcast_data")) {
                geom_col(
                  data = naive_nowcast_data,
                  aes(x = date, y = nowcast_all_naive),
                  fill = "grey", alpha = 0.5, position = "identity"
                )
              }
            } +
            {
              if (exists("naive_nowcast_known_data")) {
                geom_col(
                  data = naive_nowcast_known_data,
                  aes(x = date, y = nowcast_known_naive),
                  fill = "darkgrey", alpha = 0.5, position = "identity"
                )
              }
            } +
            geom_rect(
              data = data.frame(
                nowcast_date = forcats::fct_inorder(names(feather_days), ordered = T),
                x = current_nowcast_dates
              ),
              aes(
                xmax = x - lags[1], xmin = x - lags[length(lags)],
                ymin = minimum_y, ymax = maximum_y
              ),
              fill = "#afafb6", alpha = 0.4
            ) +
            geom_ribbon(
              aes(
                x = date,
                ymin = pmax(nowcast_all.lower, minimum_y),
                ymax = pmin(nowcast_all.upper, maximum_y),
                group = group,
                fill = model
              ),
              linetype = "dotted", alpha = 0.2
            ) +
            {
              if (exists("truth_data_cases")) {
                geom_line(
                  data = truth_data_cases,
                  aes(x = date, y = pmin(nowcast_all_true, maximum_y)),
                  color = ifelse(is.null(truth_data_index), NA, "black"),
                  alpha = 0.6, linetype = ground_truth_linetype
                )
              }
            } +
            geom_line(
              aes(
                x = date,
                y = ifelse(nowcast_all > maximum_y & lag(nowcast_all > maximum_y),
                           NA, pmin(nowcast_all, maximum_y)
                ),
                group = group, color = model
              ),
              linewidth = 0.4
            ) +
            ylab(expression(Number ~ of ~ symptom ~ onsets ~ ~ N[t])) +
            scale_x_date(
              expand = c(0, 0),
              date_labels = date_labels,
              breaks = ~ max(.x) - seq(
                7,
                feather_maxdelay + lags_list[[length(lags_list)]][1],
                by = 14
              ),
              date_breaks = date_breaks
            ) +
            scale_y_continuous(expand = expansion(add = c(0, 0))) +
            coord_cartesian(
              ylim = c(minimum_y, maximum_y), clip = "off"
            ) +
            facet_wrap(~nowcast_date, nrow = 1, scales = "free")
        }
    } else if (plot_type == "R") {
      minimum_y <- 0
      maximum_y <- 2.6
      
      PlotNowcast <- nowcast_data_cases %>%
        {
          ggplot(data = .) +
            geom_hline(yintercept = 1, linetype = "dashed", color = "grey") +
            geom_rect(
              data = data.frame(
                nowcast_date = forcats::fct_inorder(names(feather_days), ordered = T),
                x = current_nowcast_dates
              ),
              aes(
                xmax = x - lags[1],
                xmin = x - lags[length(lags)],
                ymin = 0, ymax = maximum_y
              ),
              fill = "#afafb6", alpha = 0.4
            ) +
            geom_ribbon(
              aes(
                x = date,
                ymin = pmax(R.lower, minimum_y),
                ymax = pmin(R.upper, maximum_y),
                group = group, fill = model
              ),
              linetype = "dotted", alpha = 0.2
            ) +
            {
              if (exists("truth_data_cases")) {
                geom_line(
                  data = truth_data_cases,
                  aes(x = date, y = pmin(R_true, maximum_y)),
                  color = ifelse(is.null(truth_data_index), NA, "black"),
                  alpha = 0.6, linetype = ground_truth_linetype
                )
              }
            } +
            geom_line(
              aes(
                x = date,
                y = ifelse(R > maximum_y, NA, R),
                group = group, color = model
              ),
              linewidth = 0.4
            ) +
            ylab(expression(Effective ~ reproduction ~ number ~ ~ R[t])) +
            scale_x_date(
              expand = c(0, 0),
              date_labels = date_labels,
              breaks = ~ max(.x) - seq(
                7,
                feather_maxdelay +
                  lags_list[[length(lags_list)]][1],
                by = 14
              ),
              date_breaks = date_breaks
            ) +
            scale_y_continuous(expand = expansion(add = c(0, 0))) +
            coord_cartesian(
              ylim = c(minimum_y, maximum_y), clip = "off"
            ) +
            facet_wrap(~nowcast_date, nrow = 1, scales = "free")
        }
    } else {
      stop("Unknown plot type.")
    }
    
    nowcast_guide <- guide_legend(direction = "vertical")
    
    PlotNowcast <- PlotNowcast +
      # this vertical line marks the nowcast date
      geom_vline(
        data = data.frame(
          nowcast_date = forcats::fct_inorder(names(feather_days), ordered = T),
          x = current_nowcast_dates
        ),
        aes(xintercept = x), linetype = "dotted", color = "darkgrey"
      ) +
      scale_color_manual(
        name = categories_name,
        values = model_colors[
          names(model_colors) %in% model_names[models_select]
        ]
      ) +
      scale_fill_manual(
        name = categories_name,
        values = model_colors[
          names(model_colors) %in% model_names[models_select]
        ]
      ) +
      theme_bw() +
      theme(
        legend.position = "top",
        legend.background = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_rect(
          fill = NA, color = "black", linewidth = 1
        ),
      ) +
      ggtitle(names(lags_list)[[i]]) +
      guides(color = nowcast_guide, fill = nowcast_guide)
    
    plotlegend <- cowplot::get_legend(PlotNowcast)
    PlotNowcast <- PlotNowcast +
      theme(
        legend.position = "none",
        plot.margin = unit(c(2, 2, 2, 2), "mm")
      )
    
    # plot performance
    if (is.null(performance_metric)) {
      intermediate_plots[[i]] <- cowplot::plot_grid(
        PlotNowcast,
        ncol = 1, align = "v"
      )
    } else {
      if (performance_metric == "crps") {
        plotPerformance <- plot_performance_select(m_res_list,
                                                   ifelse(plot_type == "Cases", "nowcast_all", "R"),
                                                   models_select,
                                                   ifelse(plot_type == "Cases", "mean_crps_scaled", "mean_wis"),
                                                   ifelse(plot_type == "Cases",
                                                          expression(bar(CRPS[scaled])(hat(N)[t])),
                                                          expression(bar(CRPS)(hat(R)[t]))
                                                   ),
                                                   delays_select = lags
        ) +
          coord_cartesian(xlim = c(-0.4, 0.4), expand = FALSE)
      } else if (performance_metric == "wis") {
        plotPerformance <- plot_performance_wis_select(
          m_res_list,
          ifelse(plot_type == "Cases", "nowcast_all", "R"),
          models_select,
          delays_select = lags
        )
      } else {
        stop(paste("Metric", performance_metric, "is not yet supported."))
      }
      
      plotPerformance <- plotPerformance + theme(
        strip.text = element_blank(),
        strip.background = element_blank(),
        plot.margin = unit(c(0, 2, 0.3, 2), "mm")
      )
      
      if (performance_metric == "wis" && percent_best) {
        plotPerformanceBest <- plot_performance_wis_best(
          m_res_list,
          ifelse(plot_type == "Cases", "nowcast_all", "R"),
          models_select,
          delays_select = lags
        ) +
          theme(
            strip.text = element_blank(),
            strip.background = element_blank(),
            plot.margin = unit(c(0, 2, 4, 2), "mm")
          )
      }
      
      
      # combine plots
      g <- ggplotGrob(
        ggplot() +
          theme_void() +
          annotate("polygon",
                   x = c(0, 0.317, 0.5, 1),
                   y = c(0, 1, 1, 0),
                   fill = "#afafb6",
                   alpha = 0.4
          )
      )
      PlotNowcast <- PlotNowcast +
        annotation_custom(
          g,
          xmin = -Inf,
          xmax = Inf,
          ymin = ifelse(plot_type == "Cases",
                        -(maximum_y - minimum_y) / 6.5, -0.48
          ),
          ymax = 0
        )
      
      if (performance_metric == "wis" && percent_best) {
        intermediate_plots[[i]] <- cowplot::plot_grid(
          PlotNowcast, plotPerformance, plotPerformanceBest,
          ncol = 1, rel_heights = c(0.73, 0.22, 0.09), align = "v"
        )
      } else {
        intermediate_plots[[i]] <- cowplot::plot_grid(
          PlotNowcast, plotPerformance,
          ncol = 1, rel_heights = c(0.75, 0.25), align = "v"
        )
      }
    }
  }
  
  if (is.null(performance_metric)) {
    final_plot <- cowplot::plot_grid(
      cowplot::plot_grid(NULL,
                         cowplot::plot_grid(onset_legend,
                                            plotlegend,
                                            nrow = 1,
                                            rel_widths = legends.rel_widths,
                                            align = "h"
                         ),
                         NULL,
                         nrow = 1, rel_widths = c(0.1, 0.8, 0.1)
      ),
      cowplot::plot_grid(
        plotlist = intermediate_plots,
        ncol = 1, labels = "AUTO"
      ),
      nrow = 2, rel_heights = legend_plot.rel_heights
    )
  } else if (performance_metric == "wis") {
    final_plot <- cowplot::plot_grid(
      cowplot::plot_grid(NULL,
                         cowplot::plot_grid(onset_legend,
                                            plotlegend,
                                            pattern_legend,
                                            nrow = 1,
                                            rel_widths = legends.rel_widths,
                                            align = "h"
                         ),
                         NULL,
                         nrow = 1, rel_widths = c(0.1, 0.8, 0.1)
      ),
      cowplot::plot_grid(
        plotlist = intermediate_plots,
        ncol = 1, labels = "AUTO", align = "v"
      ),
      nrow = 2, rel_heights = legend_plot.rel_heights
    )
  } else {
    final_plot <- cowplot::plot_grid(
      cowplot::plot_grid(NULL,
                         cowplot::plot_grid(onset_legend,
                                            plotlegend,
                                            nrow = 1,
                                            rel_widths = legends.rel_widths,
                                            align = "h"
                         ),
                         NULL,
                         nrow = 1, rel_widths = c(0.1, 0.8, 0.1)
      ),
      cowplot::plot_grid(
        plotlist = intermediate_plots,
        ncol = 1, labels = "AUTO"
      ),
      nrow = 2, rel_heights = legend_plot.rel_heights
    )
  }
  return(final_plot)
}
