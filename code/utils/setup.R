###########################################################################
###########################################################################
###                                                                     ###
###                                SETUP                                ###
###                                                                     ###
###########################################################################
###########################################################################

## ---------------------------------------------------------------
##                        Library imports                       -
## ---------------------------------------------------------------

here::i_am("code/utils/setup.R")

if (!exists("CLUSTER_MODE")) CLUSTER_MODE <- FALSE

if (!CLUSTER_MODE) library(qs)
if (!CLUSTER_MODE) library(nleqslv)

# library(surveillance)

# data manipulation
library(dplyr)
library(tidyr)
if (!CLUSTER_MODE) library(readr)
if (!CLUSTER_MODE) library(purrr)
if (!CLUSTER_MODE) library(tibble)
library(magrittr)
library(lubridate)
library(stringr)
if (!CLUSTER_MODE) library(forcats)
if (!CLUSTER_MODE) library(zoo)

# stats
library(extraDistr)

# stan libraries
library(cmdstanr)
library(tidybayes)

if (!CLUSTER_MODE) library(scoringutils)

# epi
library(EpiEstim)

# plotting
if (!CLUSTER_MODE) library(ggplot2)
if (!CLUSTER_MODE) library(ggdist)
if (!CLUSTER_MODE) library(ggpattern)
if (!CLUSTER_MODE) library(ggridges)
if (!CLUSTER_MODE) library(plotly)

## ---------------------------------------------------------------
##                          Constants                           -
## ---------------------------------------------------------------

Sys.setlocale("LC_ALL", "en_US")

if (!CLUSTER_MODE) {
  local_config <- here::here("code", "utils", "local_config.R")
  if (file.exists(local_config)) source(local_config)
}

additional_data <- here::here("data", "additional_data_CHE.R")
if (file.exists(additional_data)) source(additional_data)

## ---------------------------------------------------------------
##                            Utils                             -
## ---------------------------------------------------------------

source(here::here("code", "utils", "utils.R"))
source(here::here("code", "utils", "cmdstan_model_optional_profiling.R"))
source(here::here("code", "utils", "utils_summarize.R"))
if (!CLUSTER_MODE) source(here::here("code", "utils", "utils_plot.R"))
source(here::here("code", "utils", "utils_validate.R"))


source(here::here("code", "utils", "utils_standata.R"))
source(here::here("code", "utils", "utils_model.R"))
source(here::here("code", "utils", "utils_priors.R"))
source(here::here("code", "utils", "utils_fit.R"))
source(here::here("code", "utils", "utils_job.R"))
source(here::here("code", "utils", "utils_nowcast.R"))
