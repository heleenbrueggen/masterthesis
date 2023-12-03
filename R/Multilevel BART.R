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
stan4bart_fit <- stan4bart(y ~ bart(x4 + x5 + x6 + x7 + z1 + z2) + (1 + x1 + x2 + x3|group) + x1 + x2 + x3,
                           data = simdatasets$simdata_ngroup_30_groupsize_5_icc_0.5_mar_mcar_0_g_0.2[[1]])
fitted(stan4bart_fit, type = 'ranef', sample = 'train')
fitted(stan4bart_fit, type = 'indiv.bart', sample = 'train')
fitted(stan4bart_fit, type = 'indiv.ranef', sample = 'train')
fitted(stan4bart_fit, type = 'ev', sample = 'train')

# plotting y against predicted y from stan4bart
pred <- tibble(y = simdatasets$simdata_ngroup_30_groupsize_5_icc_0.5_mar_mcar_0_g_0.2[[1]]$y,
               ypred = fitted(stan4bart_fit, type = 'ev', sample = 'train'))

ggplot(data = pred, mapping = aes(
  x = y,
  y = ypred
)) +
  geom_point(color = 'royalblue', size = .5) +
  theme_minimal() +
  geom_abline(slope = 1, intercept = 0, color = 'royalblue4')

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
dbarts_fit <- rbart_vi(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + z1 + z2,
                       group.by = group,
                       data = simdatasets$simdata_ngroup_30_groupsize_5_icc_0.5_mar_mcar_0_g_0.2[[1]],
                       keepTrees = TRUE,
                       verbose = FALSE)

# Convergence
plot(dbarts_fit)

# Extracting trees
trees <- dbarts::extract(dbarts_fit, type = 'trees', sample = 'train')
fitted(dbarts_fit, type = 'bart', sample = 'train')
fitted(dbarts_fit, type = 'ranef', sample = 'train')


# plotting y against predicted y from dbarts
pred2 <- tibble(y = simdatasets$simdata_ngroup_30_groupsize_5_icc_0.5_mar_mcar_0_g_0.2[[1]]$y,
               ypred = fitted(dbarts_fit, type = 'ev', sample = 'train'))

ggplot(data = pred2, mapping = aes(
  x = y,
  y = ypred
)) +
  geom_point(color = 'royalblue', size = .5) +
  theme_minimal() +
  geom_abline(slope = 1, intercept = 0, color = 'royalblue4')

# dbarts for all data sets
dbarts_fits <- list()
for (i in 1:length(simdatasets)) {
  dbarts_fit <- simdatasets[[i]] %>%
    map(~.x %$%
          rbart_vi(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + z1 + z2,
                   group.by = group,
                   keepTrees = TRUE,
                   verbose = FALSE) %>%
          fitted(., type = 'ev', sample = 'train')) %>%
    unlist() %>%
    matrix(., ncol = length(simdatasets[[i]]))

  name <- paste('dbarts_ev', colnames(combinations)[1], combinations[i,1], colnames(combinations)[2], combinations[i,2], colnames(combinations)[3], combinations[i,3], colnames(combinations)[4], combinations[i,4], colnames(combinations)[5], combinations[i,5], sep = '_')
  dbarts_fits[[name]] <- dbarts_fit
}

dbarts_ranefs <- list()
for (i in 1:length(simdatasets)) {
  dbarts_ranef <- simdatasets[[i]] %>%
    map(~.x %$%
          rbart_vi(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + z1 + z2,
                   group.by = group,
                   keepTrees = TRUE,
                   verbose = FALSE) %>%
          fitted(., type = 'ranef', sample = 'train')) %>%
    unlist() %>%
    matrix(., ncol = length(simdatasets[[i]]))

  name <- paste('dbarts_ranef', colnames(combinations)[1], combinations[i,1], colnames(combinations)[2], combinations[i,2], colnames(combinations)[3], combinations[i,3], colnames(combinations)[4], combinations[i,4], colnames(combinations)[5], combinations[i,5], sep = '_')
  dbarts_ranefs[[name]] <- dbarts_ranef
}
################
# BART package #
################
BART_fit <- wbart(x.train = simdatasets$simdata_ngroup_30_groupsize_5_icc_0.5_mar_mcar_0_g_0.2[[1]][,3:11] %>% as.matrix(), y.train = simdatasets$simdata_ngroup_30_groupsize_5_icc_0.5_mar_mcar_0_g_0.2[[1]]$y)
BART_fit$yhat.train.mean # train data fits
BART_fit$treedraws$trees

BART_trees <- extractTreeData(BART_fit, simdatasets$simdata_ngroup_30_groupsize_5_icc_0.5_mar_mcar_0_g_0.2[[1]])
plotTree(treeData = BART_trees, treeNo = 1, iter = 1, plotType = "dendrogram")

# plotting y against predicted y from BART
pred3 <- tibble(y = simdatasets$simdata_ngroup_30_groupsize_5_icc_0.5_mar_mcar_0_g_0.2[[1]]$y,
               ypred = BART_fit$yhat.train.mean)

ggplot(data = pred3, mapping = aes(
  x = y,
  y = ypred
)) +
  geom_point(color = 'royalblue', size = .5) +
  theme_minimal() +
  geom_abline(slope = 1, intercept = 0, color = 'royalblue4')

# BART for all data sets
BART_fits <- list()
for (i in 1:length(simdatasets)) {
  BART_fits[i] <- simdatasets[[i]] %>%
    map(~.x %$%
          wbart(x.train = .[,3:11] %>% as.matrix(), y.train = .$y) %>%
          .$yhat.train.mean) %>%
    unlist() %>%
    matrix(., ncol = length(simdatasets[[i]]))
}
###################
# hebart packange #
###################
hebart_fit <- hebart(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + z1 + z2,
                     data = simdatasets$simdata_ngroup_30_groupsize_5_icc_0.5_mar_mcar_0_g_0.2[[1]],
                     group_variable = 'group')
diagnostics(hebart_fit)
predict_hebart(simdatasets$simdata_ngroup_30_groupsize_5_icc_0.5_mar_mcar_0_g_0.2[[1]], simdatasets$simdata_ngroup_30_groupsize_5_icc_0.5_mar_mcar_0_g_0.2[[1]]$group, hebart_fit, type = 'mean')

pred4 <- tibble(y = simdatasets$simdata_ngroup_30_groupsize_5_icc_0.5_mar_mcar_0_g_0.2[[1]]$y,
                ypred = hebart_fit[["y_hat"]][1, ])

ggplot(data = pred4, mapping = aes(
  x = y,
  y = ypred
)) +
  geom_point(color = 'royalblue', size = .5) +
  theme_minimal() +
  geom_abline(slope = 1, intercept = 0, color = 'royalblue4')
