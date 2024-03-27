################
# PMM analyses #
################
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
#######################
# Defining parameters #
#######################
# ngroups <- c(30, 50)
# groupsizes <- c(15, 35, 50)
# iccs <- c(.2, .5)
# mar_mcar <- c("mar", "mcar")
# miss <- c(25, 50)
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
mar_mcar <- c("mar")
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
###################
# Loading results #
###################
imp_pmm <- list()
for (i in seq_len(nrow(combinations))) {
  imp_pmm[[i]] <- map(read_rds(paste("/Volumes/Heleen\ 480GB/Master\ thesis/results/imputed/pmm/results_pmm_", names[i], ".rds", sep = ""))[1:100], ~.x$imp)
}
############################
# Plan parallel processing #
############################
cl <- makeForkCluster(5)
#####################
# Obtaining results #
#####################
lmer.model <- function(x) {
    model <- with(x, lme4::lmer(y ~ 1 + x1 + x2 + x3 + x4 + x5 + x6 + x7 + z1 + z2 + x1 * z1 + x2 * z1 + x3 * z2 + (1 + x1 + x2 + x3 | group), REML = TRUE, control = lme4:::lmerControl(optimizer = "bobyqa")))
    result <- mitml::testEstimates(mice::as.mitml.result(model), extra.pars = TRUE)

    return(result)
}
analyses_pmm <- list()
for (i in seq_len(nrow(combinations))) {
    # Logging iteration
    cat("Processing iteration:", i, "\n")
    # Perform analyses
    analyses_pmm <- pblapply(imp_pmm[[i]], lmer.model, cl = cl)
    # Save results
    write_rds(analyses_pmm, paste("/Volumes/Heleen\ 480GB/Master\ thesis/results/imputed/pmm/analyses_pmm_", names[i], ".rds", sep = ""))
}
############################
# Stop parallel processing #
############################
stopCluster(cl)