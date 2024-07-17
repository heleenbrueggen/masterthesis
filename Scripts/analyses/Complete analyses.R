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
################
# Setting path # 
################
path <- "/Volumes/Heleen 480GB/MBART-MICE files/"
#######################
# Defining parameters #
#######################
combinations <- read_rds(paste(path, "data/combinations.rds", sep = ""))
#############################
# Storing names of datasets #
#############################
names <- read_rds(paste(path, "data/names.rds", sep = ""))
#############
# Load data #
#############
# Extracting imputation objects from results
simdatasets_nomiss <- list()
for (i in seq_len(nrow(combinations))) {
  simdatasets_nomiss[[i]] <- read_rds(paste(path, "data/nomissing/simdata_", names[i], ".rds", sep = ""))
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
  write_rds(analyses_nomiss, file = paste(path, "results/nomissing/analyses_nomiss_", names[i], ".rds", sep = ""))
}
############################
# Stop parallel processing #
############################
stopCluster(cl)
