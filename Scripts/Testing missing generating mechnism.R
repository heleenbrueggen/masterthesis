########################################
# Testing missing generating mechanism #
########################################
#############
# Libraries #
#############
library(tidyverse)
library(broom.mixed)
library(naniar)
library(mice)
library(pbapply)
library(parallel)
################
# Setting seed #
################
set.seed(123)
#######################
# Defining parameters #
#######################
ngroups <- c(30, 50)
groupsizes <- c(15, 50)
iccs <- c(.5)
mar_mcar <- c("mar", "mcar")
miss <- c(50)
g <- c(.5)
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
simdatasets_miss <- list()
for (i in seq_len(nrow(combinations))) {
    simdatasets_miss[[i]] <- read_rds(paste("/Volumes/Heleen\ 480GB/Master\ thesis/data/missing/simdata_miss_", names[i], ".rds", sep = "")) %>%
        map(., ~ .x %>% select(x1, x2, x3, x4, x5, x6, x7, z1, z2, y))
}
############################
# Plan parallel processing #
############################
cl <- makeForkCluster(5)
################
# Testing MCAR #
################
missingness <- list()
for (i in seq_len(nrow(combinations))) {
    missingness[[i]] <- pblapply(simdatasets_miss[[i]], mcar_test, cl = cl) %>%
        Reduce("+", .) / length(simdatasets_miss[[i]])
}
missingness <- missingness %>% map(., ~ as_tibble(.x)) # convert to tibble
missingness %>% list_rbind() # bind all tibbles together
# Save results
write_rds(missingness, "/Volumes/Heleen 480GB/Master thesis/results/diagnostics/missingness_mcar.rds")
############################
# Stop parallel processing #
############################
stopCluster(cl)