###################
# Imputation BART #
###################
#############
# Libraries #
#############
library(devtools)
# install_github("heleenbrueggen/mice@impute.mbart", force = TRUE)
library(mice)
library(furrr)
library(magrittr)
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
simdatasets_miss <- list()
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
  simdatasets_miss[[i]] <- read_rds(paste("data/missing/", name, ".rds", sep = ""))
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
##############
# Imputation #
##############
########
# bart #
########
imputed_bart <- list()
for (i in seq_len(nrow(combinations))) {
    # Logging iteration
    cat("Processing iteration:", i, "\n")
    if (combinations[i, "miss"] != 0) {
        imputed_bart[[i]] <- simdatasets_miss[[i]] %>%
            future_map(function(x) {
                mice(x,
                    method = "bart",
                    m = 5,
                    maxit = 10,
                    seed = 123
                ) |> complete()
            }, .options = furrr_options(seed = 123))
    } else {
        imputed_bart[[i]] <- simdatasets_miss[[i]]
    }
}
#############
# stan4bart #
#############
imputed_stan4bart <- list()
for (i in seq_len(nrow(combinations))) {
    # Logging iteration
    cat("Processing iteration:", i, "\n")
    if (combinations[i, "miss"] != 0) {
        imputed_stan4bart[[i]] <- simdatasets_miss[[i]] %>%
            future_map(function(x) {
                mice(x,
                    method = "2l.bart",
                    m = 5,
                    maxit = 10,
                    seed = 123
                ) |> complete()
            }, .options = furrr_options(seed = 123))
    } else {
        imputed_stan4bart[[i]] <- simdatasets_miss[[i]]
    }
}
#########
# rbart #
#########
imputed_rbart <- list()
for (i in seq_len(nrow(combinations))) {
    # Logging iteration
    cat("Processing iteration:", i, "\n")
    if (combinations[i, "miss"] != 0) {
        imputed_rbart[[i]] <- simdatasets_miss[[i]] %>%
            future_map(function(x) {
                pred <- make.predictorMatrix(x)
                pred[, "group"] <- -2
                pred["group", "group"] <- 0
                pred[, "id"] <- 0
                mice(x,
                    method = "2l.rbart",
                    pred = pred,
                    m = 5,
                    maxit = 10,
                    seed = 123
                ) |> complete()
            }, .options = furrr_options(seed = 123))
    } else {
        imputed_rbart[[i]] <- simdatasets_miss[[i]]
    }
}
