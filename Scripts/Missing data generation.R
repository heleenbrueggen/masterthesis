###########################
# Missing data generation #
###########################
#############
# Libraries #
#############
library(furrr)
library(parallel)
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
###########################
# Defining design factors #
###########################
ngroups <- c(30, 50) # Number of groups 
groupsizes <- c(15, 50) # Group sizes 
iccs <- c(.5) # Intraclass correlation coefficient
mar_mcar <- c("mar", "mcar") # Missing data mechanism
miss <- c(50) # Percentage of missing data
g <- c(.5) # Within-group effect size
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
############################
# Plan parallel processing #
############################
cores <- detectCores() - 1 # Use all cores except one
plan(multisession, workers = cores) 
##########################
# Model based simulation #
##########################
# Defining patterns for missing mechanism
patterns <- expand.grid(c(0, 1), c(0, 1), c(0, 1), c(0, 1), c(0, 1), c(0, 1), c(0, 1), c(0, 1), c(0, 1), c(0, 1)) %>%
  filter(rowSums(.) == 9 | rowSums(.) == 8 | rowSums(.) == 7 | rowSums(.) == 6 | rowSums(.) == 5) %>%
  as.matrix()
colnames(patterns) <- c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "z1", "z2", "y")
patterns <- patterns[sample(nrow(patterns), 10, replace = FALSE), ] # Sample only 10 patterns for all possibilities

# Determining the frequency of each pattern
freq <- ampute.default.freq(patterns)

# Determining the weights of each pattern
weights <- ampute.default.weights(patterns, "MAR")
colnames(weights) <- c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "z1", "z2", "y")
# Increasing weights for z2 and x4
weights[,"z2"] <- weights[,"z2"] * 1.5
weights[,"x4"] <- weights[,"x4"] * 2

# Generating missing data
for (i in seq_len(nrow(combinations))) { # For each combination ...
  # Logging iteration
  cat("Processing iteration:", i, "\n")

  # Generating missing data
  if (combinations[i, "mar_mcar"] == "mcar") { # If missing data mechanism is MCAR ...
    simdata_miss <-
      # Read complete data
      read_rds(paste("data/complete/simdata_", names[i], ".rds", sep = "")) %>%
      # Generate missing data
      future_map(function(x) {
        others <- x %>% select(-x1, -x2, -x3, -x4, -x5, -x6, -x7, -z1, -z2, -y)
        x %>%
          select(x1, x2, x3, x4, x5, x6, x7, z1, z2, y) %>%
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
      }, .options = furrr_options(seed = 123), .progress = TRUE)
  } else { # If missing data mechanism is MAR ...
    simdata_miss <-
      # Read complete data
      read_rds(paste("data/complete/simdata_", names[i], ".rds", sep = "")) %>%
      # Generate missing data
      future_map(function(x) {
        others <- x %>% select(-x1, -x2, -x3, -x4, -x5, -x6, -x7, -z1, -z2, -y)
        x %>%
          select(x1, x2, x3, x4, x5, x6, x7, z1, z2, y) %>%
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
      }, .options = furrr_options(seed = 123), .progress = TRUE)
  }
  
  # Saving data in appropriate data folder
  write_rds(simdata_miss, file = paste("data/missing/simdata_miss_", names[i], ".rds", sep = ""))
}
############################
# Stop parallel processing #
############################
plan(sequential)