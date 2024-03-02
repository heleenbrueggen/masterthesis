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
simdatasets_miss <- list()
for (i in seq_len(nrow(combinations))) {
  # Logging iteration
  cat("Processing iteration:", i, "\n")
  if (combinations[i, "miss"] != 0) {
    if (combinations[i, "mar_mcar"] == "mcar") {
      simdatasets_miss[[i]] <-
        simdatasets[[i]] %>%
        future_map(function(x) {
          id <- x$id
          group <- x$group
          x %>%
            select(-id, -group) %>%
            group_by(z1) %>%
            ampute(
              prop = (combinations[i, "miss"] * .01),
              mech = "MCAR",
              patterns = patterns,
              freq = freq
              # bycases = FALSE
            ) %>%
            .$amp %>%
            ungroup() %>%
            mutate(id = id, group = group)
        }, .options = furrr_options(seed = 123))
    } else {
      simdatasets_miss[[i]] <-
        simdatasets[[i]] %>%
        future_map(function(x) {
          id <- x$id
          group <- x$group
          x %>%
            select(-id, -group) %>%
            group_by(z1) %>%
            ampute(
              prop = (combinations[i, "miss"] * .01),
              mech = "MAR",
              type = "RIGHT",
              patterns = patterns,
              freq = freq,
              weights = weights
              # bycases = FALSE
            ) %>%
            .$amp %>% 
            ungroup() %>%
            mutate(id = id, group = group)
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
