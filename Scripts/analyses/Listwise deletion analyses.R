##############################
# Listwise deletion analyses # 
##############################
#############
# Libraries #
#############
library(lme4)
library(lmerTest)
library(tidyverse)
library(magrittr)
library(purrr)
library(dplyr)
library(pbapply)
library(parallel)
library(readr)
library(jtools)
library(broom.mixed)
################
# Setting seed # 
################
set.seed(123)
#######################
# Defining parameters #
#######################
# ngroups <- c(30, 50)
# groupsizes <- c(15, 35, 50)
# iccs <- c(.2, .5)
# mar_mcar <- c("mcar")
# miss <- c(0)
# g <- c(.2, .5)
# combinations <- expand.grid(
#   ngroup = ngroups,
#   groupsize = groupsizes,
#   icc = iccs,
#   mar_mcar = mar_mcar,
#   miss = miss,
#   g = g
# )
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
simdata_miss <- list()
for (i in seq_len(nrow(combinations))) {
  simdata_miss[[i]] <- read_rds(paste("/Volumes/Heleen\ 480GB/Master\ thesis/data/missing/simdata_miss_", names[i], ".rds", sep = ""))[1:100] %>% map(~.x %>% na.omit())
}
############################
# Plan parallel processing #
############################
cl <- makeForkCluster(5)
#######################
# Multilevel analysis # 
#######################
lmer.model.ld <- function(x) {
  model <- x %>% lme4::lmer(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + z1 + z2 + x1 * z1 + x2 * z1 + x3 * z2 + (1 + x1 + x2 + x3 | group),
    REML = TRUE,
    control = lmerControl(optimizer = "bobyqa"), data = .)
  results <- broom.mixed::tidy(model, conf.int = TRUE)

  return(results)
}
analyses_ld <- list()
for (i in seq_len(nrow(combinations))) {
  # Logging iteration
  cat("Processing iteration:", i, "\n")
  analyses_ld <- pblapply(simdata_miss[[i]], lmer.model.ld, cl = cl)
  # Saving results
  write_rds(analyses_ld, file = paste("/Volumes/Heleen\ 480GB/Master\ thesis/results/listwise/analyses_ld_", names[i], ".rds", sep = ""))
}
############################
# Stop parallel processing #
############################
stopCluster(cl)
