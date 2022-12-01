#' Create a new CmdStanModel object
#'
#' @description \if{html}{\figure{logo.png}{options: width="25px"}}
#'   Create a new [`CmdStanModel`] object from a file containing a Stan program
#'   or from an existing Stan executable. The [`CmdStanModel`] object stores the
#'   path to a Stan program and compiled executable (once created), and provides
#'   methods for fitting the model using Stan's algorithms.
#'
#'   See the `compile` and `...` arguments for control over whether and how
#'   compilation happens.
#'
#' @export
#' @param stan_file (string) The path to a `.stan` file containing a Stan
#'   program. The helper function [write_stan_file()] is provided for cases when
#'   it is more convenient to specify the Stan program as a string. If
#'   `stan_file` is not specified then `exe_file` must be specified.
#' @param exe_file (string) The path to an existing Stan model executable. Can
#'   be provided instead of or in addition to `stan_file` (if `stan_file` is
#'   omitted some `CmdStanModel` methods like `$code()` and `$print()` will not
#'   work). This argument can only be used with CmdStan 2.27+.
#' @param compile (logical) Do compilation? The default is `TRUE`. If `FALSE`
#'   compilation can be done later via the [`$compile()`][model-method-compile]
#'   method.
#' @param profile (logical) Should the model include profile statements?
#' The default is `TRUE`. If `FALSE`, all profile statements are removed from
#' the model before compilation. This may slightly improve runtime.
#' @param ... Optionally, additional arguments to pass to the
#'   [`$compile()`][model-method-compile] method if `compile=TRUE`. These
#'   options include specifying the directory for saving the executable, turning
#'   on pedantic mode, specifying include paths, configuring C++ options, and
#'   more. See [`$compile()`][model-method-compile] for details.
#'
#' @return A [`CmdStanModel`] object.
#'
#' @seealso [install_cmdstan()], [`$compile()`][model-method-compile],
#'   [`$check_syntax()`][model-method-check_syntax]
#'
#'
#' @template seealso-docs
cmdstan_model_optional_profiling <- function(stan_file = NULL, exe_file = NULL,
                                             compile = TRUE, profile = TRUE, ...) {
  if (!profile) {
    return(
      cmdstan_model_no_profiling(
        stan_file = stan_file, exe_file = exe_file, compile = compile, ...
      )
    )
  } else {
    return(
      cmdstanr::cmdstan_model(
        stan_file = stan_file, exe_file = exe_file, compile = compile, ...
      )
    )
  }
}

#' Create a new CmdStanModel object with profiling turned off
#'
#' @description \if{html}{\figure{logo.png}{options: width="25px"}}
#'   Create a new [`CmdStanModel`] object without profiling from a file
#'   containing a Stan program. A modified version of the Stan program with
#'   all profile statements removed is created and stored temporarily.
#'   The [`CmdStanModel`] object stores the path to the Stan program and
#'   the compiled executable (once created), and provides methods for fitting
#'   the model using Stan's algorithms.
#'
#'   See the `compile` and `...` arguments for control over whether and how
#'   compilation happens.
#'
#' @export
#' @param stan_file (string) The path to a `.stan` file containing a Stan
#'   program. The helper function [write_stan_file()] is provided for cases when
#'   it is more convenient to specify the Stan program as a string. Must be
#'   specified.
#' @param exe_file (string) The path to an existing Stan model executable. Can
#'   be provided in addition to `stan_file`. This argument can only be used with
#'   CmdStan 2.27+.
#' @param compile (logical) Do compilation? The default is `TRUE`. If `FALSE`
#'   compilation can be done later via the [`$compile()`][model-method-compile]
#'   method.
#' @param ... Optionally, additional arguments to pass to the
#'   [`$compile()`][model-method-compile] method if `compile=TRUE`. These
#'   options include specifying the directory for saving the executable, turning
#'   on pedantic mode, specifying include paths, configuring C++ options, and
#'   more. See [`$compile()`][model-method-compile] for details.
#'
#' @return A [`CmdStanModel`] object without profiling.
#'
#' @seealso [install_cmdstan()], [`$compile()`][model-method-compile],
#'   [`$check_syntax()`][model-method-check_syntax]
#'   #'
#'
#' @template seealso-docs
cmdstan_model_no_profiling <- function(stan_file = NULL, exe_file = NULL,
                                       compile = TRUE, ...) {
  if (is.null(stan_file)) {
    stop("Unable to remove profiling. No 'stan_file' supplied.", call. = FALSE)
  }
  args <- list(...)
  if (!any(args$include_paths == dirname(stan_file))) {
    args$include_paths <- c(args$include_paths, dirname(stan_file))
  }
  # write files with profiling removed to temp directory
  files_no_profiling <- write_stan_files_no_profiling(
    stan_file,
    args$include_paths
  )
  args$include_paths <- NULL
  # call cmdstan_model method with paths to unprofiled stan files
  return(do.call(
    cmdstanr::cmdstan_model,
    c(list(
      stan_file = files_no_profiling$model,
      exe_file = exe_file,
      include_paths = files_no_profiling$include_paths,
      compile = compile
    ), args)
  ))
}

#' Remove profile statements from a character vector representing stan code.
#'
#' @param s Character vector representing stan code.
#'
#' @return A `character` vector of the stan code without profile statements.
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
#' with all profile statements removed.
#'
#' @param stan_file The path to a .stan file containing a Stan program.
#'
#' @param include_paths (character vector) Paths to directories where Stan
#' should look for files specified in #include directives in the Stan program.
#'
#' @param target_dir The path to a directory in which the manipulated .stan
#' files without profile statements should be stored. To avoid overriding of
#' the original .stan files, this should be different from the directory of the
#' original model and the `include_paths`.
#'
#' @return A `list` containing the path to the .stan file without profiling
#' statements and the include_paths for the included .stan files without
#' profile statements.
write_stan_files_no_profiling <- function(stan_file, include_paths = NULL,
                                          target_dir = tempdir()) {
  if (target_dir == dirname(stan_file)) {
    stop("Target directory for stan file without profiling must not be ",
      "identical to directory of original stan file.",
      call. = FALSE
    )
  }
  if (!is.null(include_paths) && any(include_paths == target_dir)) {
    stop("Target directory for stan file without profiling must not be ",
      "identical to any of the include paths.",
      call. = FALSE
    )
  }

  # remove profiling from main .stan file
  code_main_model <- paste(readLines(stan_file, warn = FALSE), collapse = "\n")
  code_main_model_no_profiling <- remove_profiling(code_main_model)
  if (!dir.exists(target_dir)) {
    dir.create(target_dir, recursive = T)
  }
  main_model <- cmdstanr::write_stan_file(
    code_main_model_no_profiling,
    dir = target_dir,
    basename = basename(stan_file)
  )

  # remove profiling from included .stan files
  if (!is.null(include_paths)) {
    include_paths_no_profiling <- rep(NA, length(include_paths))
    for (i in 1:length(include_paths)) {
      include_paths_no_profiling[i] <- file.path(
        target_dir, paste0("include_", i), basename(include_paths[i])
      )
      include_files <- list.files(
        include_paths[i],
        pattern = "*.stan", recursive = TRUE
      )
      for (f in include_files) {
        include_paths_no_profiling_fdir <- file.path(
          include_paths_no_profiling[i], dirname(f)
        )
        code_include <- paste(
          readLines(file.path(include_paths[i], f), warn = FALSE),
          collapse = "\n"
        )
        code_include_paths_no_profiling <- remove_profiling(code_include)
        if (!dir.exists(include_paths_no_profiling_fdir)) {
          dir.create(include_paths_no_profiling_fdir, recursive = T)
        }
        cmdstanr::write_stan_file(
          code_include_paths_no_profiling,
          dir = include_paths_no_profiling_fdir,
          basename = basename(f)
        )
      }
    }
  } else {
    include_paths_no_profiling <- NULL
  }

  return(list(model = main_model, include_paths = include_paths_no_profiling))
}
