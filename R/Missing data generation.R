# Missing data generation

#############
# Libraries #
#############
library(furrr)
library(mice)
library(magrittr)
library(purrr)
library(dplyr)
library(tibble)
########
# Data #
########
################
# Setting seed #
################
set.seed(123)
##########################
# Model based simulation #
##########################
simdatasets_miss <- list()
for (i in seq_len(nrow(combinations))) {
  if (combinations[i, "mar_mcar"] == "mcar") {
    simdatasets_miss[[i]] <-
      simdatasets[[i]] %>%
      future_map(function(x) {
        x %>%
          ampute(
            prop = combinations[i, "miss"],
            mech = "MCAR"
          ) %>%
          .$amp
      }, .options = furrr_options(seed = 123))
  } else {
    simdatasets_miss[[i]] <-
      simdatasets[[i]] %>%
      future_map(function(x) {
        x %>%
          ampute(
            prop = combinations[i, "miss"],
            mech = "MAR", type = "RIGHT"
          ) %>%
          .$amp
      }, .options = furrr_options(seed = 123))
  }
}
####################
# MCAR missingness #
####################
mbased_MCAR <-
  simdata %>%
  furrr::future_map(function(x) {
    x %>%
      ampute(prop = .5,
             mech = "MCAR") %>% .$amp
  }, .options = furrr_options(seed = 123))
###################
# MAR missingness #
###################
mbased_MAR <-
  simdata %>%
  furrr::future_map(function(x) {
    x %>%
      ampute(prop = .5,
             mech = "MAR", type = "RIGHT") %>% .$amp
  }, .options = furrr_options(seed = 123))
