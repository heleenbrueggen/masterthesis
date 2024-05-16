#####################
# Complete analyses # 
#####################
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
###########################
# Defining design factors #
###########################
ngroups <- c(30, 50) # Number of groups 
groupsizes <- c(15, 50) # Group sizes 
iccs <- c(.5) # Intraclass correlation coefficient
mar_mcar <- c("mar", "mcar") # Missing data mechanism
miss <- c(0) # Percentage of missing data
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
#############
# Load data #
#############
# Extracting imputation objects from results
simdatasets_nomiss <- list()
for (i in seq_len(nrow(combinations))) {
  simdatasets_nomiss[[i]] <- read_rds(paste("data/nomissing/simdata_", names[i], ".rds", sep = ""))[1:100]
}
############################
# Plan parallel processing #
############################
cl <- makeForkCluster(5)
#######################
# Multilevel analysis # 
#######################
# Define model
lmer.model.nomiss <- function(x) {
  model <- x %>% lme4::lmer(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + z1 + z2 + x1 * z1 + x2 * z1 + x3 * z2 + (1 + x1 + x2 + x3 | group),
    REML = TRUE,
    control = lmerControl(optimizer = "bobyqa"), data = .)
  results <- broom.mixed::tidy(model, conf.int = TRUE)

  return(results)
}
# Perform analyses
analyses_nomiss <- list()
for (i in seq_len(nrow(combinations))) {
  # Logging iteration
  cat("Processing iteration:", i, "\n")
  analyses_nomiss <- pblapply(simdatasets_nomiss[[i]], lmer.model.nomiss, cl = cl)
  # Saving results
  write_rds(analyses_nomiss, file = paste("results/nomissing/analyses_nomiss_", names[i], ".rds", sep = ""))
}
############################
# Stop parallel processing #
############################
stopCluster(cl)
