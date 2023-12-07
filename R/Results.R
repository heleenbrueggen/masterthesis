###########
# Results #
###########
#############
# Libraries #
#############
library(magrittr)
################
# Setting seed #
################
set.seed(123)
########
# Data #
########
load("~/Desktop/Universiteit Utrecht/Methodology & Statistics for the Behavioural, Biomedical and Social Sciences/Year 2/Master Thesis/data/simdatasets.RData")
##########
# Models #
##########
load('results/stan4bart_fits.RData')
load('results/rbart_evs.RData')
load('results/rbart_ranefs.RData')
load('results/bart_fits.RData')
load('results/gbart_fits.RData')
###########
# Results #
###########
#############
# Functions #
#############
# bias
bias <- function(x, data, com = combinations) {
  bias <- list()
  for (i in 1:length(data)) {
    y <- data[[i]] %>%
      map(~.x$y) %>%
      unlist() %>%
      matrix(., ncol = length(data[[i]]))

    a <- x[[i]] - y

    name <- paste(deparse(substitute(x)), 'bias', colnames(com)[1], com[i,1], colnames(com)[2], com[i,2], colnames(com)[3], com[i,3], colnames(com)[4], com[i,4], colnames(com)[5], com[i,5], sep = '_')
    bias[[name]] <- a
  }

  avgbias <- matrix(NA, nrow = 1, ncol = length(data))
  for(i in 1:length(bias)) {
    avgbias[,i] <- bias[[i]] %>%
      colMeans() %>%
      mean()
  }
  colnames(avgbias) <- names(data)

  return(list(bias=bias, avgbias=avgbias))
}
# mse
mse <- function(x, data, com = combinations) {
  mse <- list()
  for (i in 1:length(data)) {
    y <- data[[i]] %>%
      map(~.x$y) %>%
      unlist() %>%
      matrix(., ncol = length(data[[i]]))

    a <- (x[[i]] - y)^2

    name <- paste(deparse(substitute(x)), 'mse', colnames(com)[1], com[i,1], colnames(com)[2], com[i,2], colnames(com)[3], com[i,3], colnames(com)[4], com[i,4], colnames(com)[5], com[i,5], sep = '_')
    mse[[name]] <- a
  }

  avgmse <- matrix(NA, nrow = 1, ncol = length(data))
  for(i in 1:length(mse)) {
    avgmse[,i] <- mse[[i]] %>%
      colMeans() %>%
      mean()
  }
  colnames(avgmse) <- names(data)

  return(list(mse=mse, avgmse=avgmse))
}
#############
# stan4bart #
#############
# bias
stan4bart_bias <- bias(stan4bart_fits, simdatasets)
# mse
stan4bart_mse <- mse(stan4bart_fits, simdatasets)
# Coverage
#########
# rbart #
#########
# bias
rbart_bias <- bias(rbart_evs, simdatasets)
# mse
rbart_mse <- mse(rbart_evs, simdatasets)
# Coverage
########
# bart #
########
# bias
bart_bias <- bias(bart_fits, simdatasets)
# mse
bart_mse <- mse(bart_fits, simdatasets)
# Coverage
#########
# gbart #
#########
# bias
gbart_bias <- bias(gbart_fits, simdatasets)
# mse
gbart_mse <- mse(gbart_fits, simdatasets)
# Coverage
