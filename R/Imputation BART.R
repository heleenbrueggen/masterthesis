###################
# Imputation BART #
###################
#############
# Libraries #
#############
library(devtools)
# install_github("heleenbrueggen/mice@impute.mbart")
library(mice)
library(miceadds)
library(furrr)
library(magrittr)
library(readr)
library(lme4)
library(broom.mixed)
library(tidyverse)
################
# Setting seed #
################
set.seed(123)
# Sys.setenv("OMP_THREAD_LIMIT" = 46)
#######################
# Defining parameters #
#######################
ngroups <- c(30, 50)
groupsizes <- c(15, 35, 50)
iccs <- c(.2, .5)
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
# simdatasets_miss <- list()
# for (i in seq_len(nrow(combinations))) {
#   simdatasets_miss[[i]] <- read_rds(paste("data/missing/simdata_miss_", names[i], ".rds", sep = ""))
# }
############################
# Plan parallel processing #
############################
# plan(multisession, workers = 46)
##############
# Imputation #
##############
#######
# pmm #
#######
imputed_pmm <- list()
for (i in seq_len(nrow(combinations))) {
  # Logging iteration
  cat("Processing iteration:", i, "\n")
  # Parallel processing 
  plan(multisession, workers = 35)
  # Loading data 
  simdatasets_miss <- read_rds(paste("data/missing/simdata_miss_", names[i], ".rds", sep = ""))
  # Imputed analysis
  imputed_pmm <- simdatasets_miss %>%
    future_map(function(x) {
      # Select relevant variables
      # others <- x %>% select(-group, -x1, -x2, -x3, -x4, -x5, -x6, -x7, -z1, -z2, -y)
      x <- x %>% select(group, x1, x2, x3, x4, x5, x6, x7, z1, z2, y)
      # Create predictor matrix
      pred <- make.predictorMatrix(x)
      pred[, "group"] <- 0
      # Imputation
      imp <- mice(x,
                  method = "pmm",
                  pred = pred,
                  m = 5,
                  maxit = 10,
                  print = FALSE
      )
      # Fit model
      fit <-  with(imp, lmer(y ~ 1 + x1 + x2 + x3 + x4 + x5 + x6 + x7 + z1 + z2 + x1 * z1 + x2 * z1 + x3 * z2 + (1 + x1 + x2 + x3 | group), REML = TRUE, control = lmerControl(optimizer = "bobyqa")))
      # Obtain results 
      results <- broom.mixed::tidy((pool(fit)), conf.int = TRUE)
      
      return(list(results = results, imp = imp))
    }, .options = furrr_options(seed = 123), .progress = TRUE)
  
  # Saving imputed results
  write_rds(imputed_pmm, file = paste("results/imputed/pmm/results_pmm_", names[i], ".rds", sep = ""))
  
  # Stop parallel processing 
  plan(sequential)
}
##################
# multilevel pmm #
##################
imputed_2l.pmm <- list()
for (i in seq_len(nrow(combinations))) {
  # Logging iteration
  cat("Processing iteration:", i, "\n")
  # Parallel processing 
  plan(multisession, workers = 35)
  # Loading data 
  simdatasets_miss <- read_rds(paste("data/missing/simdata_miss_", names[i], ".rds", sep = ""))
  # Imputed analysis
  imputed_2l.pmm <- simdatasets_miss %>%
    future_map(function(x) {
      # Select relevant variables
      # others <- x %>% select(-group, -x1, -x2, -x3, -x4, -x5, -x6, -x7, -z1, -z2, -y)
      x <- x %>% select(group, x1, x2, x3, x4, x5, x6, x7, z1, z2, y)
      # Add interaction variables
      x <- data.frame(x, 
                      x1.z1 = NA, x2.z1 = NA, x3.z2 = NA)
      # Create predictor matrix
      pred <- make.predictorMatrix(x)
      pred["group", ] <- 0
      pred["x1",] <- c(-2, 0, 4, 4, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1)
      pred["x2",] <- c(-2, 4, 0, 4, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1)
      pred["x3",] <- c(-2, 4, 4, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0)
      pred["x4",] <- c(-2, 4, 4, 4, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1)
      pred["x5",] <- c(-2, 4, 4, 4, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1)
      pred["x6",] <- c(-2, 4, 4, 4, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1)
      pred["x7",] <- c(-2, 4, 4, 4, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1)
      pred["z1",] <- c(-2, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1)
      pred["z2",] <- c(-2, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0)
      pred["y",] <- c(-2, 4, 4, 4, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1)
      pred["x1.z1",] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
      pred["x2.z1",] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
      pred["x3.z2",] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
      # Create methods
      meth <- make.method(x)
      meth[c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "z1", "z2", "y", "x1.z1", "x2.z1", "x3.z2")] <- c("2l.pmm", "2l.pmm", "2l.pmm", "2l.pmm", "2l.pmm", "2l.pmm", "2l.pmm", "2lonly.mean", "2lonly.mean", "2l.pmm", "~I(x1 * z1)", "~I(x2 * z1)", "~I(x3 * z2)")
      imp <- mice(x,
                  method = meth,
                  pred = pred,
                  m = 5,
                  maxit = 10,
                  print = FALSE
      )
      # Fit model
      fit <-  with(imp, lmer(y ~ 1 + x1 + x2 + x3 + x4 + x5 + x6 + x7 + z1 + z2 + x1 * z1 + x2 * z1 + x3 * z2 + (1 + x1 + x2 + x3 | group), REML = TRUE, control = lmerControl(optimizer = "bobyqa"))) 
      # Obtain results
      results <- broom.mixed::tidy((pool(fit)), conf.int = TRUE)
      
      return(list(results = results, imp = imp))
    }, .options = furrr_options(seed = 123), .progress = TRUE)
  
  # Saving imputed results
  write_rds(imputed_2l.pmm, file = paste("results/imputed/2l.pmm/results_2l.pmm_", names[i], ".rds", sep = ""))
  
  # Stop parallel processing 
  plan(sequential)
}
########
# bart #
########
imputed_bart <- list()
for (i in seq_len(nrow(combinations))) {
  # Logging iteration
  cat("Processing iteration:", i, "\n")
  # Parallel processing 
  plan(multisession, workers = 35)
  # Loading data 
  simdatasets_miss <- read_rds(paste("data/missing/simdata_miss_", names[i], ".rds", sep = ""))
  # Imputed analysis
  imputed_bart <- simdatasets_miss %>%
    future_map(function(x) {
      # Select relevant variables
      # others <- x %>% select(-group, -x1, -x2, -x3, -x4, -x5, -x6, -x7, -z1, -z2, -y)
      x <- x %>% select(group, x1, x2, x3, x4, x5, x6, x7, z1, z2, y)
      # Create predictor matrix
      pred <- make.predictorMatrix(x)
      pred[, "group"] <- 0
      imp <- mice(x,
                  method = "bart",
                  pred = pred,
                  m = 5,
                  maxit = 10,
                  print = FALSE
      ) 
      fit <-  with(imp, lmer(y ~ 1 + x1 + x2 + x3 + x4 + x5 + x6 + x7 + z1 + z2 + x1 * z1 + x2 * z1 + x3 * z2 + (1 + x1 + x2 + x3 | group), REML = TRUE, control = lmerControl(optimizer = "bobyqa"))) 
      
      results <- broom.mixed::tidy((pool(fit)), conf.int = TRUE)
      
      return(list(results = results, imp = imp))
    }, .options = furrr_options(seed = 123), .progress = TRUE)
  
  # Saving imputed results
  write_rds(imputed_bart, file = paste("results/imputed/bart/results_bart_", names[i], ".rds", sep = ""))
  
  # Stop parallel processing
  plan(sequential)
}
#########
# rbart #
#########
imputed_rbart <- list()
for (i in seq_len(nrow(combinations))) {
  # Logging iteration
  cat("Processing iteration:", i, "\n")
  # Parallel processing 
  plan(multisession, workers = 35)
  # Loading data 
  simdatasets_miss <- read_rds(paste("data/missing/simdata_miss_", names[i], ".rds", sep = ""))
  # Imputed analysis
  imputed_rbart <- simdatasets_miss %>%
    future_map(function(x) {
      # Select relevant variables
      # others <- x %>% select(-group, -x1, -x2, -x3, -x4, -x5, -x6, -x7, -z1, -z2, -y)
      x <- x %>% select(group, x1, x2, x3, x4, x5, x6, x7, z1, z2, y)
      # Create predictor matrix
      pred <- make.predictorMatrix(x)
      pred[, "group"] <- -2
      pred["group", "group"] <- 0
      imp <- mice(x,
                  method = "2l.rbart",
                  pred = pred,
                  m = 5,
                  maxit = 10,
                  print = FALSE
      )
      fit <-  with(imp, lmer(y ~ 1 + x1 + x2 + x3 + x4 + x5 + x6 + x7 + z1 + z2 + x1 * z1 + x2 * z1 + x3 * z2 + (1 + x1 + x2 + x3 | group), REML = TRUE, control = lmerControl(optimizer = "bobyqa"))) 
      
      results <- broom.mixed::tidy((pool(fit)), conf.int = TRUE)
      
      return(list(results = results, imp = imp))
    }, .options = furrr_options(seed = 123), .progress = TRUE)
  
  # Saving imputed results
  write_rds(imputed_rbart, file = paste("results/imputed/rbart/results_rbart_", names[i], ".rds", sep = ""))
  
  # Stop parallel processing
  plan(sequential)
}
#############
# stan4bart #
#############
imputed_stan4bart <- list()
for (i in seq_len(nrow(combinations))) {
  # Logging iteration
  cat("Processing iteration:", i, "\n")
  # Parallel processing 
  plan(multisession, workers = 35)
  # Loading data 
  simdatasets_miss <- read_rds(paste("data/missing/simdata_miss_", names[i], ".rds", sep = ""))
  # Imputed analysis
  imputed_stan4bart <- simdatasets_miss %>%
    future_map(function(x) {
      # Select relevant variables
      # others <- x %>% select(-group, -x1, -x2, -x3, -x4, -x5, -x6, -x7, -z1, -z2, -y)
      x <- x %>% select(group, x1, x2, x3, x4, x5, x6, x7, z1, z2, y)
      # Create predictor matrix
      pred <- make.predictorMatrix(x)
      pred["x1", ] <- c(-2, 0, 2, 2, 1, 1, 1, 1, 1, 1, 1)
      pred["x2", ] <- c(-2, 2, 0, 2, 1, 1, 1, 1, 1, 1, 1)
      pred["x3", ] <- c(-2, 2, 2, 0, 1, 1, 1, 1, 1, 1, 1)
      pred["x4", ] <- c(-2, 2, 2, 2, 0, 1, 1, 1, 1, 1, 1)
      pred["x5", ] <- c(-2, 2, 2, 2, 1, 0, 1, 1, 1, 1, 1)
      pred["x6", ] <- c(-2, 2, 2, 2, 1, 1, 0, 1, 1, 1, 1)
      pred["x7", ] <- c(-2, 2, 2, 2, 1, 1, 1, 0, 1, 1, 1)
      pred["z1", ] <- c(-2, 2, 2, 2, 1, 1, 1, 1, 0, 1, 1)
      pred["z2", ] <- c(-2, 2, 2, 2, 1, 1, 1, 1, 1, 0, 1)
      pred["y", ] <- c(-2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 0)
      pred["group", ] <- 0
      imp <- mice(x,
                  method = "2l.bart",
                  pred = pred,
                  m = 5,
                  maxit = 10,
                  print = FALSE
      )
      fit <-  with(imp, lmer(y ~ 1 + x1 + x2 + x3 + x4 + x5 + x6 + x7 + z1 + z2 + x1 * z1 + x2 * z1 + x3 * z2 + (1 + x1 + x2 + x3 | group), REML = TRUE, control = lmerControl(optimizer = "bobyqa"))) 
      
      results <- broom.mixed::tidy((pool(fit)), conf.int = TRUE)
      
      return(list(results = results, imp = imp))
    }, .options = furrr_options(seed = 123), .progress = TRUE)
  
  # Saving imputed results
  write_rds(imputed_stan4bart, file = paste("results/imputed/stan4bart/results_stan4bart_", names[i], ".rds", sep = ""))
  
  # Stop paralllel processing 
  plan(sequential)
}
############################
# Stop parallel processing #
############################
# plan(sequential)
