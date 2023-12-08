###########
# Results #
###########
#############
# Libraries #
#############
library(magrittr)
library(ggplot2)
library(tibble)
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
# Combinations for current models (excluding MAR and MCAR)
combinations2 <- expand.grid(ngroup = ngroups,
                            groupsize = groupsizes,
                            icc = iccs,
                            g = g)
#############
# Functions #
#############
# bias
bias <- function(x, data, com = combinations) {
  bias <- list()
  for (i in 1:length(data)) { # Extract y value for all data sets
    y <- data[[i]] %>%
      map(~.x$y) %>%
      unlist() %>%
      matrix(., ncol = length(data[[i]]))

    a <- x[[i]] - y
    # name <- paste('bias', colnames(combinations)[1], combinations[i,1], colnames(combinations)[2], combinations[i,2], colnames(combinations)[3], combinations[i,3], colnames(combinations)[5], combinations[i,5], sep = '_')
    bias[[i]] <- a # Store bias in a list
   }
  avgbias <- matrix(NA, nrow = length(bias), ncol = 1)
  for(i in 1:length(bias)) { # Calculate average bias per data set
    avgbias[i,] <- bias[[i]] %>%
      rowMeans() %>%
      mean()
  }
  #rownames(avgbias) <- names(bias)
  return(list(bias=bias, avgbias=avgbias))
}
# mse
mse <- function(x, data, com = combinations) {
  mse <- list()
  for (i in 1:length(data)) { # Extract y value for all data sets
    y <- data[[i]] %>%
      map(~.x$y) %>%
      unlist() %>%
      matrix(., ncol = length(data[[i]]))

    a <- (x[[i]] - y)^2
    # name <- paste('mse', colnames(combinations)[1], combinations[i,1], colnames(combinations)[2], combinations[i,2], colnames(combinations)[3], combinations[i,3], colnames(combinations)[5], combinations[i,5], sep = '_')
    mse[[i]] <- a
  }

  avgmse <- matrix(NA, nrow = length(mse), ncol = 1)
  for(i in 1:length(mse)) { # Calculate average mse per data set
    avgmse[i,] <- mse[[i]] %>%
      colMeans() %>%
      mean()
  }
  # rownames(avgmse) <- names(mse)
  return(list(mse=mse, avgmse=avgmse)) # Store output in a list
}
#############
# stan4bart #
#############
# bias
stan4bart_bias <- bias(stan4bart_fits, simdatasets)
# mse
stan4bart_mse <- mse(stan4bart_fits, simdatasets)
# coverage
#########
# rbart #
#########
# bias
rbart_bias <- bias(rbart_evs, simdatasets)
# mse
rbart_mse <- mse(rbart_evs, simdatasets)
# coverage
########
# bart #
########
# bias
bart_bias <- bias(bart_fits, simdatasets)
# mse
bart_mse <- mse(bart_fits, simdatasets)
# coverage
#########
# gbart #
#########
# bias
gbart_bias <- bias(gbart_fits, simdatasets)
# mse
gbart_mse <- mse(gbart_fits, simdatasets)
# coverage
#########
# Plots #
#########
# bias for all models and all combinations
bias_models <- tibble(stan4bart = stan4bart_bias$avgbias,
       rbart = rbart_bias$avgbias,
       bart = bart_bias$avgbias,
       gbart = gbart_bias$avgbias,
       dataset = paste(colnames(combinations)[1], combinations[,1], colnames(combinations)[2], combinations[,2], colnames(combinations)[3], combinations[,3], colnames(combinations)[5], combinations[,5], sep = '_')) %>%
  cbind(., combinations2) %>%
  pivot_longer(cols = c(stan4bart, rbart, bart, gbart), names_to = 'model', values_to = 'bias') %>%
  group_by(dataset, model, icc) %>%
  summarise(across(bias, ~ mean(.x)), .groups = "drop")
# mse for all models and all combinations
mse_models <- tibble(stan4bart = stan4bart_mse$avgmse,
                      rbart = rbart_mse$avgmse,
                      bart = bart_mse$avgmse,
                      gbart = gbart_mse$avgmse,
                      dataset = paste(colnames(combinations)[1], combinations[,1], colnames(combinations)[2], combinations[,2], colnames(combinations)[3], combinations[,3], colnames(combinations)[5], combinations[,5], sep = '_')) %>%
  cbind(., combinations2) %>%
  pivot_longer(cols = c(stan4bart, rbart, bart, gbart), names_to = 'model', values_to = 'mse') %>%
  group_by(dataset, model, icc) %>%
  summarise(across(mse, ~ mean(.x)), .groups = "drop")
# color palette
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# bias plots
bias_plot <- function (biasdata, icc_value) {
  data <- biasdata %>%
    dplyr::filter(icc == icc_value)

  ggplot(data ,aes(
      x = dataset,
      y = bias,
      fill = model,
      color = model
    )) +
    geom_bar(stat = 'identity', position = 'identity', width = .05) +
    geom_point(stat = 'identity') +
    labs(
      title = paste0('Average bias for all models, icc =', icc),
      x = 'Simulated dataset',
      y = 'Bias') +
    theme_bw() +
    theme(axis.text.x = element_text(size = 8, angle = 90, hjust = 1)) +
    scale_fill_manual(values=cbbPalette) +
    scale_color_manual(values=cbbPalette)
}
bias_plot(bias_models, icc_value = 0)
bias_plot(bias_models, icc_value = 0.05)
bias_plot(bias_models, icc_value = 0.3)
bias_plot(bias_models, icc_value = 0.5)
# mse plots
mse_plot <- function (msedata, icc_value) {
  data <- msedata %>%
    dplyr::filter(icc == icc_value)

  ggplot(data ,aes(
    x = dataset,
    y = mse,
    fill = model,
    color = model
  )) +
    geom_bar(stat = 'identity', position = 'identity', width = .05) +
    geom_point(stat = 'identity') +
    labs(
      title = paste0('Average MSE for all models, icc =', icc),
      x = 'Simulated dataset',
      y = 'MSE') +
    theme_bw() +
    theme(axis.text.x = element_text(size = 8, angle = 90, hjust = 1)) +
    scale_fill_manual(values=cbbPalette) +
    scale_color_manual(values=cbbPalette)
}
mse_plot(mse_models, icc_value = 0)
mse_plot(mse_models, icc_value = 0.05)
mse_plot(mse_models, icc_value = 0.3)
mse_plot(mse_models, icc_value = 0.5)
