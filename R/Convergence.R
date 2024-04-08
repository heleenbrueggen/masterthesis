######################
# 2l.pmm convergence #
######################
#############
# Libraries #
#############
library(tidyverse)
library(broom.mixed)
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
mar_mcar <- c("mar", "mcar")
miss <- c(50)
g <- c(.5)
combinations <- expand.grid(
  ngroup = ngroups
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
imp_2l.pmm <- list()
for (i in seq_len(nrow(combinations))) {
  imp_2l.pmm[[i]] <- map(read_rds(paste("/Volumes/Heleen\ 480GB/Master\ thesis/results/imputed/2l.pmm/results_2l.pmm_", names[i], ".rds", sep = ""))[1:100], ~.x$imp)
}
imp_bart <- list()
for (i in seq_len(nrow(combinations))) {
  imp_bart[[i]] <- map(read_rds(paste("/Volumes/Heleen\ 480GB/Master\ thesis/results/imputed/bart/results_bart_", names[i], ".rds", sep = ""))[1:100], ~.x$imp)
}
imp_rbart <- list()
for (i in seq_len(nrow(combinations))) {
  imp_rbart[[i]] <- map(read_rds(paste("/Volumes/Heleen\ 480GB/Master\ thesis/results/imputed/rbart/results_rbart_", names[i], ".rds", sep = ""))[1:100], ~.x$imp)
}
imp_stan4bart <- list()
for (i in seq_len(nrow(combinations))) {
  imp_stan4bart[[i]] <- map(read_rds(paste("/Volumes/Heleen\ 480GB/Master\ thesis/results/imputed/stan4bart/results_stan4bart_", names[i], ".rds", sep = ""))[1:20], ~.x$imp)
}
############################
# Plan parallel processing #
############################
cl <- makeForkCluster(5)
#####################
# Convergence plots #
#####################
plot.imp <- function(imp, layout) {
    plot <- mice:::plot.mids(imp, layout = c(4, 7))
    stripplot <- mice:::stripplot.mids(imp, layout = c(4, 7))
    return(list(convergence = plot, stripplot = stripplot))
}
# pmm 
plots_pmm <- list()
for (i in seq_len(nrow(combinations))) {
    # Logging iteration
    cat("Processing iteration:", i, "\n")
    # Plots 
    plots_pmm[[i]] <- pblapply(imp_pmm[[i]], plot.imp, cl = cl)
}
# 2l.pmm
plots_2l.pmm <- list()
for (i in seq_len(nrow(combinations))) {
    # Logging iteration
    cat("Processing iteration:", i, "\n")
    # Plots 
    plots_2l.pmm[[i]] <- pblapply(imp_2l.pmm[[i]], plot.imp, cl = cl)
}
# bart
plots_bart <- list()
for (i in seq_len(nrow(combinations))) {
    # Logging iteration
    cat("Processing iteration:", i, "\n")
    # Plots 
    plots_bart[[i]] <- pblapply(imp_bart[[i]], plot.imp, cl = cl)
}
# rbart 
plots_rbart <- list()
for (i in seq_len(nrow(combinations))) {
    # Logging iteration
    cat("Processing iteration:", i, "\n")
    # Plots 
    plots_rbart[[i]] <- pblapply(imp_rbart[[i]], plot.imp, cl = cl)
}
# stan4bart
plots_stan4bart <- list()
for (i in seq_len(nrow(combinations))) {
    # Logging iteration
    cat("Processing iteration:", i, "\n")
    # Plots 
    plots_stan4bart[[i]] <- pblapply(imp_stan4bart[[i]], plot.imp, cl = cl)
}
############################
# Stop parallel processing #
############################
stopCluster(cl)