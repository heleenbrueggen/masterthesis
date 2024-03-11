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
groupsizes <- c(15, 35, 50)
iccs <- c(0, .3)
mar_mcar <- c("mar", "mcar")
miss <- c(25, 50)
g <- c(.2, .5)
combinations <- expand.grid(
  ngroup = ngroups,
  groupsize = groupsizes,
  icc = iccs,
  mar_mcar = mar_mcar,
  miss = miss,
  g = g
)
#############################
# Storing names of datasets #
#############################
names <- rep(NA, nrow(combinations))
for (i in seq_len(nrow(combinations))) {
  names[i] <- paste(
    colnames(combinations)[1], combinations[i, 1],
    colnames(combinations)[2], combinations[i, 2],
    colnames(combinations)[3], combinations[i, 3],
    colnames(combinations)[4], combinations[i, 4],
    colnames(combinations)[5], combinations[i, 5],
    colnames(combinations)[6], combinations[i, 6],
    sep = "_"
  )
}
#############
# Load data #
#############
simdatasets <- list()
for (i in seq_len(nrow(combinations))) {
  simdatasets[[i]] <- read_rds(paste("data/complete/simdata_", names[i], ".rds", sep = ""))
}
##########################
# Model based simulation #
##########################
# Defining patterns for missing mechanism
patterns <- expand.grid(c(0, 1), c(0, 1), c(0, 1), c(0, 1), c(0, 1), c(0, 1), c(0, 1), c(0, 1), c(0, 1), c(0, 1)) %>%
  filter(rowSums(.) == 9 | rowSums(.) == 8 | rowSums(.) == 7 | rowSums(.) == 6 | rowSums(.) == 5) %>%
  as.matrix()
colnames(patterns) <- c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "z1", "z2", "y")

# Determining the frequency of each pattern
freq <- ampute.default.freq(patterns)
# freq[rowSums(patterns) == 9] <- 0.4 / sum(rowSums(patterns) == 9)
# freq[rowSums(patterns) == 8] <- 0.2 / sum(rowSums(patterns) == 8)
# freq[rowSums(patterns) == 7] <- 0.2 / sum(rowSums(patterns) == 7)
# freq[rowSums(patterns) == 6] <- 0.1 / sum(rowSums(patterns) == 6)
# freq[rowSums(patterns) == 5] <- 0.1 / sum(rowSums(patterns) == 5)

# Determining the weights of each pattern
weights <- ampute.default.weights(patterns, "MAR")
colnames(weights) <- c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "z1", "z2", "y")
weights[,"z2"] <- weights[,"z2"] * 1.5
weights[,"x4"] <- weights[,"x4"] * 2

# test <- ampute(simdatasets[[129]][[1]][,3:12], prop = .03333, patterns = patterns, freq = freq, weights = weights, mech = "MAR", bycases = FALSE)$amp

# Generating missing data
for (i in seq_len(nrow(combinations))) {
  # Logging iteration
  cat("Processing iteration:", i, "\n")
  if (combinations[i, "mar_mcar"] == "mcar") {
    simdata_miss <-
      simdatasets[[i]] %>%
      future_map(function(x) {
        others <- x %>% select(-x1, -x2, -x3, -x4, -x5, -x6, -x7, -z1, -z2, -y)
        x %>%
          select(x1, x2, x3, x4, x5, x6, z7, z1, z2, y) %>%
          group_by(z1) %>%
          ampute(
            prop = (combinations[i, "miss"] * .01),
            mech = "MCAR",
            patterns = patterns,
            freq = freq
          ) %>%
          .$amp %>%
          ungroup() %>%
          cbind(., others)
      }, .options = furrr_options(seed = 123))
  } else {
    simdata_miss <-
      simdatasets[[i]] %>%
      future_map(function(x) {
        others <- x %>% select(-x1, -x2, -x3, -x4, -x5, -x6, -x7, -z1, -z2, -y)
        x %>%
          select(x1, x2, x3, x4, x5, x6, z7, z1, z2, y) %>%
          group_by(z1) %>%
          ampute(
            prop = (combinations[i, "miss"] * .01),
            mech = "MAR",
            type = "RIGHT",
            patterns = patterns,
            freq = freq,
            weights = weights
          ) %>%
          .$amp %>%
          ungroup() %>%
          cbind(., others)
      }, .options = furrr_options(seed = 123))
  }

  # Saving data in appropriate data folder
  write_rds(simdata_miss, file = paste("data/missing/simdata_miss_", names[i], ".rds", sep = ""))
}
