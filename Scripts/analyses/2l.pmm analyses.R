###################
# 2l.pmm analyses #
###################
#############
# Libraries #
#############
library(tidyverse)
library(broom.mixed)
library(lme4)
library(mitml)
library(mice)
library(pbapply)
library(parallel)
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
###################
# Loading results #
###################
# Extracting imputation objects from results 
imp_2l.pmm <- list()
for (i in seq_len(nrow(combinations))) {
  imp_2l.pmm[[i]] <- map(read_rds(paste(path, "results/imputed/2l.pmm/results_2l.pmm_", names[i], ".rds", sep = "")), ~.x$imp)
}
############################
# Plan parallel processing #
############################
cl <- makeForkCluster(5)
#####################
# Obtaining results #
#####################
# Define model
lmer.model <- function(x) {
    model <- with(x, lme4::lmer(y ~ 1 + x1 + x2 + x3 + x4 + x5 + x6 + x7 + z1 + z2 + x1 * z1 + x2 * z1 + x3 * z2 + (1 + x1 + x2 + x3 | group), REML = TRUE, control = lme4:::lmerControl(optimizer = "bobyqa")))
    result <- mitml::testEstimates(mice::as.mitml.result(model), extra.pars = TRUE)
    return(result)
}
# Perform analyses
analyses_2l.pmm <- list()
for (i in seq_len(nrow(combinations))) {
    # Logging iteration
    cat("Processing iteration:", i, "\n")
    # Perform analyses
    analyses_2l.pmm <- pblapply(imp_2l.pmm[[i]], lmer.model, cl = cl)
    # Save results
    write_rds(analyses_2l.pmm, paste(path, "results/imputed/2l.pmm/analyses_2l.pmm_", names[i], ".rds", sep = ""))
}
############################
# Stop parallel processing #
############################
stopCluster(cl)
