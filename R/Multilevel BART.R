###################
# Multilevel BART #
###################
#############
# Libraries #
#############
library(stan4bart)
library(BART)
library(dbarts)
library(bartMachine)
library(bartMan)
library(tidyverse)
library(ggplot2)
library(magrittr)
library(hebartBase)
################
# Setting seed #
################
set.seed(123)
#####################
# stan4bart package #
#####################
# stan4bart for all data sets
stan4bart_fits <- list()
for (i in 1:length(simdatasets)) {
  stan4bart_fits[i] <- simdatasets[[i]] %>%
    map(~.x %$%
          stan4bart(y ~ bart(x4 + x5 + x6 + x7 + z1 + z2) + (1 + x1 + x2 + x3|group) + x1 + x2 + x3) %>%
          fitted(., type = 'ev', sample = 'train')) %>%
    unlist() %>%
    matrix(., ncol = length(simdatasets[[i]]))
}
##################
# dbarts package #
##################
# rbart_vi for all data sets
rbart_ev <- list()
for (i in 1:length(simdatasets)) {
  rbartr_ev <- simdatasets[[i]] %>%
    map(~.x %$%
          rbart_vi(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + z1 + z2,
                   group.by = group,
                   verbose = FALSE) %>%
          fitted(., type = 'ev', sample = 'train')) %>%
    unlist() %>%
    matrix(., ncol = length(simdatasets[[i]]))

  name <- paste('rbart_ev', colnames(combinations)[1], combinations[i,1], colnames(combinations)[2], combinations[i,2], colnames(combinations)[3], combinations[i,3], colnames(combinations)[4], combinations[i,4], colnames(combinations)[5], combinations[i,5], sep = '_')
  rbart_ev[[name]] <- rbart_ev
}

rbart_ranefs <- list()
for (i in 1:length(simdatasets)) {
  rbart_ranef <- simdatasets[[i]] %>%
    map(~.x %$%
          rbart_vi(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + z1 + z2,
                   group.by = group,
                   verbose = FALSE) %>%
          fitted(., type = 'ranef', sample = 'train')) %>%
    unlist() %>%
    matrix(., ncol = length(simdatasets[[i]]))

  name <- paste('rbart_ranef', colnames(combinations)[1], combinations[i,1], colnames(combinations)[2], combinations[i,2], colnames(combinations)[3], combinations[i,3], colnames(combinations)[4], combinations[i,4], colnames(combinations)[5], combinations[i,5], sep = '_')
  rbart_ranefs[[name]] <- rbart_ranef
}

# single level bart model
bart_fits <- list()
for (i in 1:length(simdatasets)) {
  bart_fit <- simdatasets[[i]] %>%
    map(~.x %$%
          bart2(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + z1 + z2,
                verbose = FALSE) %>%
          fitted(., type = 'ev', sample = 'train')) %>%
    unlist() %>%
    matrix(., ncol = length(simdatasets[[i]]))

  name <- paste('bart_ev', colnames(combinations)[1], combinations[i,1], colnames(combinations)[2], combinations[i,2], colnames(combinations)[3], combinations[i,3], colnames(combinations)[4], combinations[i,4], colnames(combinations)[5], combinations[i,5], sep = '_')
  bart_fits[[name]] <- bart_fit
}

# single level BART model with group dummy
gbart_fits <- list()
for (i in 1:length(simdatasets)) {
  gbart_fit <- simdatasets[[i]] %>%
    map(~.x %$%
          bart2(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + z1 + z2 + group,
                verbose = FALSE) %>%
          fitted(., type = 'ev', sample = 'train')) %>%
    unlist() %>%
    matrix(., ncol = length(simdatasets[[i]]))

  name <- paste('gbart_ev', colnames(combinations)[1], combinations[i,1], colnames(combinations)[2], combinations[i,2], colnames(combinations)[3], combinations[i,3], colnames(combinations)[4], combinations[i,4], colnames(combinations)[5], combinations[i,5], sep = '_')
  gbart_fits[[name]] <- gbart_fit
}
