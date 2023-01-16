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
library(tidyverse)
library(magrittr)
library(lubridate)
library(stringr)
library(zoo)

# stats
library(extraDistr)

# stan libraries
library(cmdstanr)
library(tidybayes)

library(scoringutils)

# epi
library(EpiEstim)

# plotting
library(plotly)
library(ggdist)
if (!CLUSTER_MODE) library(ggridges)

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
source(here::here("code", "utils", "utils_plot.R"))
source(here::here("code", "utils", "utils_validate.R"))

source(here::here("code", "make_nowcast.R"))
