###################
# Imputation BART #
###################
#############
# Libraries #
#############
library(devtools)
# Install the version of mice including necessary imputation methods
# install_github("heleenbrueggen/mice@impute.mbart")
library(mice)
library(miceadds)
library(furrr)
library(magrittr)
library(readr)
library(lme4)
library(dplyr)
library(broom.mixed)
library(tidyverse)
library(doParallel)
library(parallel)
library(pbapply)
library(stan4bart)
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
############################
# Plan parallel processing #
############################
cores <- detectCores() - 1 # Use all cores except one
cl <- makeCluster(cores)
##############
# Imputation #
##############
#######
# pmm #
#######
# Function for pmm imputation
pmm.analysis <- function(x) {
  # Select relevant variables 
  x <- x |> dplyr::select(group, x1, x2, x3, x4, x5, x6, x7, z1, z2, y)
  # Create predictor matrix
  pred <- mice::make.predictorMatrix(x)
  pred[, "group"] <- 0
  # Imputation
  imp <- mice::mice(x,
                    method = "pmm",
                    pred = pred,
                    m = 5,
                    maxit = 10,
                    print = FALSE)

  return(list(imp = imp))
}
# Imputation pmm
for (i in seq_len(nrow(combinations))) { # For each combination ...
  # Logging iteration
  cat("Processing iteration:", i, "\n")
  # Loading data
  simdatasets_miss <- read_rds(paste(path, "data/missing/simdata_miss_", names[i], ".rds", sep = ""))[1:100] %>%
    map(., \(x) x$data) # Only use 100 datasets
  # Imputation
  imputed_pmm <- pblapply(simdatasets_miss, pmm.analysis, cl = cl)

  # Saving imputed results
  write_rds(imputed_pmm, file = paste(path, "results/imputed/pmm/results_pmm_", names[i], ".rds", sep = ""))
}
##################
# multilevel pmm #
##################
# Function for multilevel pmm imputation
pmm2l.analysis <- function(x) {
  # Select relevant variables
  # others <- x %>% select(-group, -x1, -x2, -x3, -x4, -x5, -x6, -x7, -z1, -z2, -y)
  x <- x |> dplyr::select(group, x1, x2, x3, x4, x5, x6, x7, z1, z2, y)
  # Add interaction variables
  x <- data.frame(x, 
                  x1.z1 = NA, x2.z1 = NA, x3.z2 = NA)
  # Create predictor matrix
  pred <- mice::make.predictorMatrix(x)
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
  library(miceadds)
  meth <- mice::make.method(x)
  meth[c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "z1", "z2", "y", "x1.z1", "x2.z1", "x3.z2")] <- c("2l.pmm", "2l.pmm", "2l.pmm", "2l.pmm", "2l.pmm", "2l.pmm", "2l.pmm", "2lonly.mean", "2lonly.mean", "2l.pmm", "~I(x1 * z1)", "~I(x2 * z1)", "~I(x3 * z2)")
  # Imputation
  imp <- mice::mice(x,
                    method = meth,
                    pred = pred,
                    m = 5,
                    maxit = 10,
                    print = FALSE
  )
  
  return(list(imp = imp))
}
# Imputation 2l.pmm
for (i in seq_len(nrow(combinations))) { # For each combination ...
  # Logging iteration
  cat("Processing iteration:", i, "\n")
  # Loading data 
  simdatasets_miss <- read_rds(paste(path, "data/missing/simdata_miss_", names[i], ".rds", sep = ""))[1:100] %>%
    map(., \(x) x$data) # Only use 100 datasets
  # Imputation
  imputed_2l.pmm <- pblapply(simdatasets_miss, pmm2l.analysis, cl = cl)
  
  # Saving imputed results
  write_rds(imputed_2l.pmm, file = paste(path, "results/imputed/2l.pmm/results_2l.pmm_", names[i], ".rds", sep = ""))
}
########
# bart #
########
# Function for bart imputation
bart.analysis <- function(x) {
  # Select relevant variables
  x <- x |> dplyr::select(group, x1, x2, x3, x4, x5, x6, x7, z1, z2, y)
  # Create predictor matrix
  pred <- mice::make.predictorMatrix(x)
  pred[, "group"] <- 0
  # Imputation
  imp <- mice::mice(x,
                    method = "bart",
                    pred = pred,
                    m = 5,
                    maxit = 10,
                    print = FALSE
  )
  
  return(list(imp = imp))
}
for (i in seq_len(nrow(combinations))) { # For each combination ...
  # Logging iteration
  cat("Processing iteration:", i, "\n")
  # Loading data 
  simdatasets_miss <- read_rds(paste(path, "data/missing/simdata_miss_", names[i], ".rds", sep = ""))[1:100] %>%
    map(., \(x) x$data) # Only use 100 datasets
  # Imputation
  imputed_bart <- pblapply(simdatasets_miss, bart.analysis, cl = cl)
  
  # Saving imputed results
  write_rds(imputed_bart, file = paste(path, "results/imputed/bart/results_bart_", names[i], ".rds", sep = ""))
}
#########
# rbart #
#########
# Function for rbart imputation
rbart.analysis <- function(x) {
  # Select relevant variables
  x <- x |> dplyr::select(group, x1, x2, x3, x4, x5, x6, x7, z1, z2, y)
  # Create predictor matrix
  pred <- mice::make.predictorMatrix(x)
  pred[, "group"] <- -2
  pred["group", "group"] <- 0
  # Imputation
  imp <- mice::mice(x,
                    method = "2l.rbart",
                    pred = pred,
                    m = 5,
                    maxit = 10,
                    print = FALSE
  )
  
  return(list(imp = imp))
}
for (i in seq_len(nrow(combinations))) { # For each combination ...
  # Logging iteration
  cat("Processing iteration:", i, "\n")
  # Loading data 
  simdatasets_miss <- read_rds(paste(path, "data/missing/simdata_miss_", names[i], ".rds", sep = ""))[1:100] %>%
    map(., \(x) x$data) # Only use 100 datasets
  # Imputation
  imputed_rbart <- pblapply(simdatasets_miss, rbart.analysis, cl = cl)
  
  # Saving imputed results
  write_rds(imputed_rbart, file = paste(path, "results/imputed/rbart/results_rbart_", names[i], ".rds", sep = ""))
}
#############
# stan4bart #
#############
# Function for stan4bart imputation
stan4bart.analysis <- function(x) {
  # Select relevant variables
  x <- x |> dplyr::select(group, x1, x2, x3, x4, x5, x6, x7, z1, z2, y)
  # Create predictor matrix
  pred <- mice::make.predictorMatrix(x)
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
  # Imputation
  imp <- mice::mice(x,
                    method = "2l.bart",
                    pred = pred,
                    m = 5,
                    maxit = 10,
                    print = FALSE
  )
  
  return(list(imp = imp))
}
for (i in seq_len(nrow(combinations))) { # For each combination ...
  # Logging iteration
  cat("Processing iteration:", i, "\n")
  # Loading data 
  simdatasets_miss <- read_rds(paste(path, "data/missing/simdata_miss_", names[i], ".rds", sep = ""))[1:100] %>%
    map(., \(x) x$data) # Only use 100 datasets
  # Imputed analysis
  imputed_stan4bart <- pblapply(simdatasets_miss, stan4bart.analysis, cl = cl)
  
  # Saving imputed results
  write_rds(imputed_stan4bart, file = paste(path, "results/imputed/stan4bart/results_stan4bart_", names[i], ".rds", sep = ""))
}
############################
# Stop parallel processing #
############################
stopCluster(cl)