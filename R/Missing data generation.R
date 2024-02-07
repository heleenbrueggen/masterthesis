###########################
# Missing data generation #
###########################
#############
# Libraries #
#############
library(furrr)
library(mice)
library(magrittr)
library(purrr)
library(dplyr)
library(tibble)
library(readr)
################
# Setting seed #
################
set.seed(123)
#######################
# Defining parameters #
#######################
ngroups <- c(30, 50)
groupsizes <- c(5, 15, 35, 50)
iccs <- c(0, .05, .3, .5)
mar_mcar <- c("mar", "mcar")
miss <- c(0, 25, 50)
g <- c(.2, .5, .8)
combinations <- expand.grid(
  ngroup = ngroups,
  groupsize = groupsizes,
  icc = iccs,
  mar_mcar = mar_mcar,
  miss = miss,
  g = g
)
#############
# Load data #
#############
simdatasets <- list()
for (i in seq_len(nrow(combinations))) {
  name <- paste("simdata",
    colnames(combinations)[1], combinations[i, 1],
    colnames(combinations)[2], combinations[i, 2],
    colnames(combinations)[3], combinations[i, 3],
    colnames(combinations)[4], combinations[i, 4],
    colnames(combinations)[5], combinations[i, 5],
    colnames(combinations)[6], combinations[i, 6],
    sep = "_"
  )
  simdatasets[[i]] <- read_rds(paste("data/complete/", name, ".rds", sep = ""))
}
#############################
# Storing names of datasets #
#############################
names <- rep(NA, 576)
for (i in seq_len(nrow(combinations))) {
  names[i] <- paste("simdata",
    colnames(combinations)[1], combinations[i, 1],
    colnames(combinations)[2], combinations[i, 2],
    colnames(combinations)[3], combinations[i, 3],
    colnames(combinations)[4], combinations[i, 4],
    colnames(combinations)[5], combinations[i, 5],
    colnames(combinations)[6], combinations[i, 6],
    sep = "_"
  )
}
##########################
# Model based simulation #
##########################
simdatasets_miss <- list()
for (i in seq_len(nrow(combinations))) {
  if (combinations[i, "miss"] != 0) {
    if (combinations[i, "mar_mcar"] == "mcar") {
      simdatasets_miss[[i]] <-
        simdatasets[[i]] %>%
        future_map(function(x) {
          id <- x$id
          group <- x$group
          x %>%
            select(-id, -group) %>%
            ampute(
              prop = combinations[i, "miss"],
              mech = "MCAR"
            ) %>%
            .$amp %>%
            mutate(
              id = id,
              group = group
            )
        }, .options = furrr_options(seed = 123))
    } else {
      simdatasets_miss[[i]] <-
        simdatasets[[i]] %>%
        future_map(function(x) {
          id <- x$id
          group <- x$group
          x %>%
            select(-id, -group) %>%
            ampute(
              prop = combinations[i, "miss"],
              mech = "MAR", type = "RIGHT"
            ) %>%
            .$amp %>%
            mutate(
              id = id,
              group = group
            )
        }, .options = furrr_options(seed = 123))
    }
  } else {
    simdatasets_miss[[i]] <- simdatasets[[i]]
  }
}
#######################
# Saving missing data # 
#######################
write_rds(simdatasets_miss, file = "data/missing/simdatasets_miss.rds")
for (i in seq_len(nrow(combinations))) {
  name <- paste("simdata_miss",
    colnames(combinations)[1], combinations[i, 1],
    colnames(combinations)[2], combinations[i, 2],
    colnames(combinations)[3], combinations[i, 3],
    colnames(combinations)[4], combinations[i, 4],
    colnames(combinations)[5], combinations[i, 5],
    colnames(combinations)[6], combinations[i, 6],
    sep = "_"
  )
  write_rds(simdatasets_miss[[i]], file = paste("data/missing/", name, ".rds", sep = ""))
}
