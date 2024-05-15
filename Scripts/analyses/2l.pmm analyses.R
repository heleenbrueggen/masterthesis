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
###########################
# Defining design factors #
###########################
ngroups <- c(30, 50) # Number of groups 
groupsizes <- c(15, 50) # Group sizes 
iccs <- c(.5) # Intraclass correlation coefficient
mar_mcar <- c("mar", "mcar") # Missing data mechanism
miss <- c(50) # Percentage of missing data
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
###################
# Loading results #
###################
# Extracting imputation objects from results 
imp_2l.pmm <- list()
for (i in seq_len(nrow(combinations))) {
  imp_2l.pmm[[i]] <- map(read_rds(paste("results/imputed/2l.pmm/results_2l.pmm_", names[i], ".rds", sep = ""))[1:100], ~.x$imp)
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
    write_rds(analyses_2l.pmm, paste("results/imputed/2l.pmm/analyses_2l.pmm_", names[i], ".rds", sep = ""))
}
############################
# Stop parallel processing #
############################
stopCluster(cl)
