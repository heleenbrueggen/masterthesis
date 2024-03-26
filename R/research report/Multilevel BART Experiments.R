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
bla <- stan4bart(y ~ bart(x4 + x5 + x6 + x7 + z1 + z2 + x1 + x2 + x3) + (1 + x1 + x2 + x3|group),
                           data = simdatasets$simdata_ngroup_30_groupsize_5_icc_0.5_mar_mcar_0_g_0.2[[1]],
                           bart_args = list(verbose = FALSE,
                                            k = 2.0,
                                            n.samples = 1500L,
                                            n.burn = 1500L,
                                            n.thin = 5L))
fitted(stan4bart_fit, type = 'ranef', sample = 'train')
fitted(stan4bart_fit, type = 'indiv.bart', sample = 'train')
fitted(stan4bart_fit, type = 'indiv.ranef', sample = 'train')
fitted(stan4bart_fit, type = 'ev', sample = 'train')

# plotting y against predicted y from stan4bart
pred <- tibble(y = simdatasets$simdata_ngroup_30_groupsize_5_icc_0.5_mar_mcar_0_g_0.2[[1]]$y,
               ypred = fitted(bla, type = 'ev', sample = 'train'))

ggplot(data = pred, mapping = aes(
  x = y,
  y = ypred
)) +
  geom_point(color = 'royalblue', size = .5) +
  theme_minimal() +
  geom_abline(slope = 1, intercept = 0, color = 'royalblue4')
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
################
# BART package #
################
BART_fit <- wbart(x.train = simdatasets$simdata_ngroup_30_groupsize_5_icc_0.5_mar_mcar_0_g_0.2[[1]][,3:11] %>% as.matrix(), y.train = simdatasets$simdata_ngroup_30_groupsize_5_icc_0.5_mar_mcar_0_g_0.2[[1]]$y)
BART_fit$yhat.train.mean # train data fits
BART_fit$treedraws$trees

BART_trees <- extractTreeData(BART_fit, simdatasets$simdata_ngroup_30_groupsize_5_icc_0.5_mar_mcar_0_g_0.2[[1]])
plotTree(treeData = BART_trees, treeNo = 1, iter = 1, plotType = "dendrogram")

# plotting y against predicted y from BART
pred34 <- tibble(y = simdatasets$simdata_ngroup_30_groupsize_5_icc_0.5_mar_mcar_0_g_0.2[[1]]$y,
                 ypred = BART_fit$yhat.train.mean)

ggplot(data = pred4, mapping = aes(
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

pred5 <- tibble(y = simdatasets$simdata_ngroup_30_groupsize_5_icc_0.5_mar_mcar_0_g_0.2[[1]]$y,
                ypred = hebart_fit[["y_hat"]][1, ])

ggplot(data = pred5, mapping = aes(
  x = y,
  y = ypred
)) +
  geom_point(color = 'royalblue', size = .5) +
  theme_minimal() +
  geom_abline(slope = 1, intercept = 0, color = 'royalblue4')

