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
################
# Setting seed #
################
set.seed(123)
#####################
# stan4bart package #
#####################
stan4bart_fit <- stan4bart(y ~ bart(z1*(x1 + x2) + x3*z2 + x4 + x5 + x6 + x7) + (x1|group) + (x2|group) + (x3|group), data = simdatasets$simdata_ngroup_30_groupsize_5_icc_0.5_mar_mcar_0_g_0.2[[1]])
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
trees <- extract(dbarts_fit, type = 'trees', sample = 'train')

# plotting y against predicted y from dbarts
pred <- tibble(y = simdatasets$simdata_ngroup_30_groupsize_5_icc_0.5_mar_mcar_0_g_0.2[[1]]$y,
               ypred = fitted(dbarts_fit, type = 'ev', sample = 'train'))

ggplot(data = pred, mapping = aes(
  x = y,
  y = ypred
)) +
  geom_point(color = 'royalblue', size = .5) +
  theme_minimal() +
  geom_abline(slope = 1, intercept = 0, color = 'royalblue4') +
  geom_rug(color = 'royalblue', linewidth = .1)

# dbarts for all data sets
dbarts_fits <- list()
for (i in 1:length(simdatasets)) {
  dbarts_fits[i] <- simdatasets[[i]] %>%
    map(~.x %$%
          rbart_vi(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + z1 + z2,
                   group.by = group,
                   keepTrees = TRUE,
                   verbose = FALSE) %>%
          fitted(., type = 'ev', sample = 'train')) %>%
    Reduce("+", .) / length(simdatasets[[i]])
}
################
# BART package #
################
BART_fit <- wbart(x.train = simdatasets$simdata_ngroup_30_groupsize_5_icc_0.5_mar_mcar_0_g_0.2[[1]][,3:11] %>% as.matrix(), y.train = simdatasets$simdata_ngroup_30_groupsize_5_icc_0.5_mar_mcar_0_g_0.2[[1]]$y)
BART_fit$yhat.train.mean # train data fits
BART_fit$treedraws$trees

BART_trees <- extractTreeData(BART_fit, simdatasets$simdata_ngroup_30_groupsize_5_icc_0.5_mar_mcar_0_g_0.2[[1]])
plotTree(treeData = BART_trees, treeNo = 1, iter = 1, plotType = "dendrogram")

# Obtaining the tree from BART model
BART_trees <- utils::read.table(text = BART_fit$treedraws$trees,
                                skip = 1,
                                fill = NA,
                                col.names = c("node", "var", "splitValue", "leafValue"))
BART_trees$var <- names(BART_fit$varcount.mean)[BART_trees$var + 1] # as vars are indexed at 0
BART_trees$splitID <- BART_trees$splitValue + 1
BART_trees$tier <- as.integer(floor(log2(BART_trees$node)))


# plotting y against predicted y from BART
pred <- tibble(y = simdatasets$simdata_ngroup_30_groupsize_5_icc_0.5_mar_mcar_0_g_0.2[[1]]$y,
               ypred = BART_fit$yhat.train.mean)

ggplot(data = pred, mapping = aes(
  x = y,
  y = ypred
)) +
  geom_point(color = 'royalblue', size = .5) +
  theme_minimal() +
  geom_abline(slope = 1, intercept = 0, color = 'royalblue4') +
  geom_rug(color = 'royalblue', linewidth = .1)

