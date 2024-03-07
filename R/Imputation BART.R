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
# pred <- make.predictorMatrix(simdatasets_miss[[65]][[1]])
# pred[, "id"] <- 0
# imp <- mice(simdatasets_miss[[66]][[1]], method = "pmm", m = 5, maxit = 10, seed = 123)
# fit <- with(imp, lmer(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + z1 + z2 + x1 * z1 + x2 * z1 + x3 * z2 + (x1 + x2 + x3| group), REML = FALSE))
# summary(pool(fit))
# testEstimates(as.mitml.result(fit), var.comp = TRUE)$var.comp
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
    if (combinations[i, "miss"] != 0) {
        imputed_pmm[[i]] <- simdatasets_miss[[i]] %>%
            future_map(function(x) {
                pred <- make.predictorMatrix(x)
                pred[, "id"] <- 0
                imp <- mice(x,
                    method = "pmm",
                    pred = pred,
                    m = 5,
                    maxit = 10,
                    seed = 123
                ) #|> complete()
                fit <-  with(imp, lmer(y ~ 1 + x1 + x2 + x3 + x4 + x5 + x6 + x7 + z1 + z2 + x1 * z1 + x2 * z1 + x3 * z2 + (1 + x1 + x2 + x3 + x4 + x5 + x6| group), REML = FALSE)) 
                summary(pool(fit))
                testEstimates(as.mitml.result(fit), var.comp = TRUE)$var.comp
            }, .options = furrr_options(seed = 123))
    } else {
        imputed_pmm[[i]] <- simdatasets_miss[[i]]
    }
}
##################
# multilevel pmm #
##################
# In progress
imputed_2l.pmm <- list()
for (i in seq_len(nrow(combinations))) {
    # Logging iteration
    cat("Processing iteration:", i, "\n")
    if (combinations[i, "miss"] != 0) {
        imputed_2l.pmm[[i]] <- simdatasets_miss[[i]] %>%
            future_map(function(x) {
                pred <- make.predictorMatrix(x)
                # pred[, "group"] <- -2
                pred[, "id"] <- 0
                pred["id", "id"] <- 0
                pred["x1",] <- c(0, -2, 0, 2, 2, 2, 2, 2, 1, 1, 1, 1)
                pred["x2",] <- c(0, -2, 2, 0, 2, 2, 2, 2, 1, 1, 1, 1)
                pred["x3",] <- c(0, -2, 2, 2, 0, 2, 2, 2, 1, 1, 1, 1)
                pred["x4",] <- c(0, -2, 2, 2, 2, 0, 2, 2, 1, 1, 1, 1)
                pred["x5",] <- c(0, -2, 2, 2, 2, 2, 0, 2, 1, 1, 1, 1)
                pred["x6",] <- c(0, -2, 2, 2, 2, 2, 2, 0, 1, 1, 1, 1)
                pred["x7",] <- c(0, -2, 2, 2, 2, 2, 2, 2, 0, 1, 1, 1)
                pred["z1",] <- c(0, -2, 2, 2, 2, 2, 2, 2, 1, 0, 1, 1)
                pred["z2",] <- c(0, -2, 2, 2, 2, 2, 2, 2, 1, 1, 0, 1)
                pred["y",] <- c(0, -2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 0)
                meth <- make.method(x)
                meth[c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "z1", "z2", "y")] <- c("2l.pmm", "2l.pmm", "2l.pmm", "2l.pmm", "2l.pmm", "2l.pmm", "2l.pmm", "2lonly.pmm", "2lonly.pmm", "2l.pmm")
                mice(x,
                    method = meth,
                    pred = pred,
                    m = 5,
                    maxit = 10,
                    seed = 123
                ) |> complete()
            }, .options = furrr_options(seed = 123))
    } else {
        imputed_2l.pmm[[i]] <- simdatasets_miss[[i]]
    }
}
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
                    maxit = 5,
                    seed = 123
                ) |> complete()
            }, .options = furrr_options(seed = 123))
    } else {
        imputed_bart[[i]] <- simdatasets_miss[[i]]
    }
}
write_rds(imputed_bart, file = "data/imputed/imputed_bart.rds")
imputed_bart <- read_rds("data/imputed/imputed_bart.rds")
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
write_rds(imputed_rbart, file = "data/imputed/imputed_rbart.rds")
imputed_rbart <- read_rds("data/imputed/imputed_rbart.rds")
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
                pred <- make.predictorMatrix(x)
                pred[, "id"] <- 0
                pred["x1", ] <- c(0, 2, 2, 2, 2, 2, 1, 1, 1, 1, 0, -2)
                pred["x2", ] <- c(2, 0, 2, 2, 2, 2, 1, 1, 1, 1, 0, -2)
                pred["x3", ] <- c(2, 2, 0, 2, 2, 2, 1, 1, 1, 1, 0, -2)
                pred["x4", ] <- c(2, 2, 2, 0, 2, 2, 1, 1, 1, 1, 0, -2)
                pred["x5", ] <- c(2, 2, 2, 2, 0, 2, 1, 1, 1, 1, 0, -2)
                pred["x6", ] <- c(2, 2, 2, 2, 2, 0, 1, 1, 1, 1, 0, -2)
                pred["x7", ] <- c(2, 2, 2, 2, 2, 2, 0, 1, 1, 1, 0, -2)
                pred["z1", ] <- c(2, 2, 2, 2, 2, 2, 1, 0, 1, 1, 0, -2)
                pred["z2", ] <- c(2, 2, 2, 2, 2, 2, 1, 1, 0, 1, 0, -2)
                pred["y", ] <- c(2, 2, 2, 2, 2, 2, 1, 1, 1, 0, 0, -2)
                pred["group", ] <- 0
                pred["id", ] <- 0
                meth <- make.method(x)
                meth[c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "z1", "z2", "y")] <- "2l.bart"
                mice(x,
                    method = meth,
                    pred = pred,
                    m = 5,
                    maxit = 10,
                    seed = 123
                ) |> complete()
            }, .options = furrr_options(seed = 123))
    } else {
        imputed_stan4bart[[i]] <- simdatasets_miss[[i]]
    }
}