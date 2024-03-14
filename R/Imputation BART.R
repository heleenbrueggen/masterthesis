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
simdatasets_miss <- list()
for (i in seq_len(nrow(combinations))) {
  simdatasets_miss[[i]] <- read_rds(paste("data/missing/simdata_miss_", names[i], ".rds", sep = ""))
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
        imputed_pmm <- simdatasets_miss[[i]] %>%
            future_map(function(x) {
                # Select relevant variables
                # others <- x %>% select(-group, -x1, -x2, -x3, -x4, -x5, -x6, -x7, -z1, -z2, -y)
                x <- x %>% select(group, x1, x2, x3, x4, x5, x6, x7, z1, z2, y) # Group eruit halen
                # Create predictor matrix
                pred <- make.predictorMatrix(x)
                pred[, "group"] <- 0
                # Imputation
                imp <- mice(x,
                    method = "pmm",
                    pred = pred,
                    m = 5,
                    maxit = 10
                )
                # Fit model
                fit <-  with(imp, lmer(y ~ 1 + x1 + x2 + x3 + x4 + x5 + x6 + x7 + z1 + z2 + x1 * z1 + x2 * z1 + x3 * z2 + (1 + x1 + x2 + x3 | group), REML = FALSE))
                # Obtain results 
                results <- broom.mixed::tidy((pool(fit)), conf.int = TRUE)
            }, .options = furrr_options(seed = 123))
    
    # Saving imputed results
    write_rds(imputed_pmm, file = paste("results/imputed/pmm/results_pmm_", names[i], ".rds", sep = ""))
}
##################
# multilevel pmm #
##################
# In progress
imputed_2l.pmm <- list() # testen, covergence plots 
for (i in seq_len(nrow(combinations))) {
    # Logging iteration
    cat("Processing iteration:", i, "\n")
        imputed_2l.pmm <- simdatasets_miss[[i]] %>%
            future_map(function(x) {
                # Select relevant variables
                # others <- x %>% select(-group, -x1, -x2, -x3, -x4, -x5, -x6, -x7, -z1, -z2, -y)
                x <- x %>% select(group, x1, x2, x3, x4, x5, x6, x7, z1, z2, y)
                # Add interaction variables
                x <- data.frame(x, 
                x1.z1 = NA, x2.z1 = NA, x3.x2 = NA)
                # Create predictor matrix
                pred <- make.predictorMatrix(x)
                pred["group", ] <- 0
                pred["x1",] <- c(-2, 0, 4, 4, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1) # Geen random slopes voor x variabelen
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
                pred["x2.z3",] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
                # Create methods
                meth <- make.method(x)
                meth[c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "z1", "z2", "y", "x1.z1", "x2.z1", "x3.z2")] <- c("2l.pmm", "2l.pmm", "2l.pmm", "2l.pmm", "2l.pmm", "2l.pmm", "2l.pmm", "2lonly.pmm", "2lonly.pmm", "2l.pmm", "~I(x1 * z1)", "~I(x2 * z1)", "~I(x3 * z2)")
                imp <- mice(x,
                    method = meth,
                    pred = pred,
                    m = 5,
                    maxit = 10
                )
                # Fit model
                fit <-  with(imp, lmer(y ~ 1 + x1 + x2 + x3 + x4 + x5 + x6 + x7 + z1 + z2 + x1 * z1 + x2 * z1 + x3 * z2 + (1 + x1 + x2 + x3 | group), REML = FALSE)) 
                # Obtain results
                broom.mixed::tidy((pool(fit)), conf.int = TRUE)
            }, .options = furrr_options(seed = 123))
    
    # Saving imputed results
    write_rds(imputed_2l.pmm, file = paste("results/imputed/2l.pmm/results_2l.pmm_", names[i], ".rds", sep = ""))
}
########
# bart #
########
imputed_bart <- list()
for (i in seq_len(nrow(combinations))) {
    # Logging iteration
    cat("Processing iteration:", i, "\n")
        imputed_bart <- simdatasets_miss[[i]] %>%
            future_map(function(x) {
                # Select relevant variables
                # others <- x %>% select(-group, -x1, -x2, -x3, -x4, -x5, -x6, -x7, -z1, -z2, -y)
                x <- x %>% select(group, x1, x2, x3, x4, x5, x6, x7, z1, z2, y)
                # Create predictor matrix
                pred <- make.predictorMatrix(x)
                pred[, "group"] <- 0
                imp <- mice(x,
                    method = "bart",
                    m = 5,
                    maxit = 10
                ) 
                fit <-  with(imp, lmer(y ~ 1 + x1 + x2 + x3 + x4 + x5 + x6 + x7 + z1 + z2 + x1 * z1 + x2 * z1 + x3 * z2 + (1 + x1 + x2 + x3 | group), REML = FALSE)) 

                broom.mixed::tidy((pool(fit)), conf.int = TRUE)
            }, .options = furrr_options(seed = 123))

    # Saving imputed results
    write_rds(imputed_bart, file = paste("results/imputed/bart/results_bart_", names[i], ".rds", sep = ""))
}
#########
# rbart #
#########
imputed_rbart <- list()
for (i in seq_len(nrow(combinations))) {
    # Logging iteration
    cat("Processing iteration:", i, "\n")
        imputed_rbart <- simdatasets_miss[[i]] %>%
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
                    maxit = 10
                )
                fit <-  with(imp, lmer(y ~ 1 + x1 + x2 + x3 + x4 + x5 + x6 + x7 + z1 + z2 + x1 * z1 + x2 * z1 + x3 * z2 + (1 + x1 + x2 + x3 | group), REML = FALSE)) 

                broom.mixed::tidy((pool(fit)), conf.int = TRUE)
            }, .options = furrr_options(seed = 123))
    
    # Saving imputed results
    write_rds(imputed_rbart, file = paste("results/imputed/rbart/results_rbart_", names[i], ".rds", sep = ""))
}
#############
# stan4bart #
#############
imputed_stan4bart <- list()
for (i in seq_len(nrow(combinations))) {
    # Logging iteration
    cat("Processing iteration:", i, "\n")
        imputed_stan4bart <- simdatasets_miss[[i]] %>%
            future_map(function(x) {
                # Select relevant variables
                # others <- x %>% select(-group, -x1, -x2, -x3, -x4, -x5, -x6, -x7, -z1, -z2, -y)
                x <- x %>% select(group, x1, x2, x3, x4, x5, x6, x7, z1, z2, y)
                # Create predictor matrix
                pred <- make.predictorMatrix(x)
                pred["x1", ] <- c(-2, 0, 2, 2, 1, 1, 1, 1, 1, 1, 1) # Geen random slopes voor x variabelen
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
                meth <- make.method(x)
                meth[c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "z1", "z2", "y")] <- "2l.bart"
                imp <- mice(x,
                    method = meth,
                    pred = pred,
                    m = 5,
                    maxit = 10
                )
                fit <-  with(imp, lmer(y ~ 1 + x1 + x2 + x3 + x4 + x5 + x6 + x7 + z1 + z2 + x1 * z1 + x2 * z1 + x3 * z2 + (1 + x1 + x2 + x3 | group), REML = FALSE)) 

                broom.mixed::tidy((pool(fit)), conf.int = TRUE)
            }, .options = furrr_options(seed = 123))
    
    # Saving imputed results
    write_rds(imputed_stan4bart, file = paste("results/imputed/stan4bart/results_stan4bart_", names[i], ".rds", sep = ""))
}