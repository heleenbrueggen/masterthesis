###################
# Multilevel BART #
###################
#############
# Libraries #
#############
library(stan4bart)
library(dbarts)
library(tidyverse)
library(magrittr)
library(doParallel)
################
# Setting seed #
################
set.seed(123)
########
# Data #
########
load("~/Desktop/Universiteit Utrecht/Methodology & Statistics for the Behavioural, Biomedical and Social Sciences/Year 2/Master Thesis/data/simdatasets.RData")
#####################
# stan4bart package #
#####################
# stan4bart for all data sets
stan4bart_fits <- list()
for (i in 1:length(simdatasets)) {
  stan4bart_fit <- simdatasets[[i]] %>%
    map(~.x %$%
          stan4bart(y ~ bart(x4 + x5 + x6 + x7 + z1 + z2) + (1 + x1 + x2 + x3|group) + x1 + x2 + x3,
                    cores = detectCores() - 1,
                    verbose = -1,
                    bart_args = list(k = 2.0,
                                     n.samples = 1500L,
                                     n.burn = 1500L,
                                     n.thin = 5L)) %>%
          fitted(., type = 'ev', sample = 'train')) %>%
    unlist() %>%
    matrix(., ncol = length(simdatasets[[i]]))

  name <- paste('stan4bart_ev', colnames(combinations)[1], combinations[i,1], colnames(combinations)[2], combinations[i,2], colnames(combinations)[3], combinations[i,3], colnames(combinations)[4], combinations[i,4], colnames(combinations)[5], combinations[i,5], sep = '_')
  stan4bart_fits[[name]] <- stan4bart_fit
}
##################
# dbarts package #
##################
# rbart_vi for all data sets
rbart_evs <- list()
rbart_ranefs <- list()
for (i in 1:length(simdatasets)) {
  rbart <- simdatasets[[i]] %>%
    map(~.x %$%
          rbart_vi(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + z1 + z2,
                   group.by = group,
                   verbose = FALSE))

  rbart_ev <- rbart %>%
    map(~.x %$%
          fitted(., type = 'ev', sample = 'train')) %>%
    unlist() %>%
    matrix(., ncol = length(simdatasets[[i]]))

  rbart_ranef <- rbart %>%
    map(~.x %$%
          fitted(., type = 'ranef', sample = 'train')) %>%
    unlist() %>%
    matrix(., ncol = length(simdatasets[[i]]))

  name_ev <- paste('rbart_ev', colnames(combinations)[1], combinations[i,1], colnames(combinations)[2], combinations[i,2], colnames(combinations)[3], combinations[i,3], colnames(combinations)[4], combinations[i,4], colnames(combinations)[5], combinations[i,5], sep = '_')
  rbart_evs[[name_ev]] <- rbart_ev

  name_ranef <- paste('rbart_ranef', colnames(combinations)[1], combinations[i,1], colnames(combinations)[2], combinations[i,2], colnames(combinations)[3], combinations[i,3], colnames(combinations)[4], combinations[i,4], colnames(combinations)[5], combinations[i,5], sep = '_')
  rbart_ranefs[[name_ranef]] <- rbart_ranef
}
# single level bart model
bart_fits <- list()
for (i in 1:length(simdatasets)) {
  bart_fit <- simdatasets[[i]] %>%
    map(~.x %$%
          bart2(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + z1 + z2,
                verbose = FALSE,
                k = 2.0,
                n.samples = 1500L,
                n.burn = 1500L,
                n.thin = 5L) %>%
          fitted(., type = 'ev', sample = 'train')) %>%
    unlist() %>%
    matrix(., ncol = length(simdatasets[[i]]))

  name <- paste('bart_ev', colnames(combinations)[1], combinations[i,1], colnames(combinations)[2], combinations[i,2], colnames(combinations)[3], combinations[i,3], colnames(combinations)[4], combinations[i,4], colnames(combinations)[5], combinations[i,5], sep = '_')
  bart_fits[[name]] <- bart_fit
}
# single level bart model with group dummy
gbart_fits <- list()
for (i in 1:length(simdatasets)) {
  gbart_fit <- simdatasets[[i]] %>%
    map(~.x %$%
          bart2(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + z1 + z2 + group,
                verbose = FALSE,
                k = 2.0,
                n.samples = 1500L,
                n.burn = 1500L,
                n.thin = 5L) %>%
          fitted(., type = 'ev', sample = 'train')) %>%
    unlist() %>%
    matrix(., ncol = length(simdatasets[[i]]))

  name <- paste('gbart_ev', colnames(combinations)[1], combinations[i,1], colnames(combinations)[2], combinations[i,2], colnames(combinations)[3], combinations[i,3], colnames(combinations)[4], combinations[i,4], colnames(combinations)[5], combinations[i,5], sep = '_')
  gbart_fits[[name]] <- gbart_fit
}
##################
# Saving results #
##################
save(stan4bart_fits, file = 'results/stan4bart_fits.RData')
save(rbart_evs, file = 'results/rbart_evs.RData')
save(rbart_ranefs, file = 'results/rbart_ranefs.RData')
save(bart_fits, file = 'results/bart_fits.RData')
save(gbart_fits, file = 'results/gbart_fits.RData')
