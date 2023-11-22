###################
# Multilevel BART #
###################
#############
# Libraries #
#############
library(stan4bart)
library(BART)
library(lme4)
library(dbarts)
library(bartMachine)
library(bartMan)
library(tidyverse)
library(ggplot2)
################
# Setting seed #
################
set.seed(123)
################
# lme4 package #
################
lme4_fit <- lmer(y ~ 1 + x1 + x2 + x3 + x4 + x5 + x6 + x7 + z1 + z2 + x1*z1 + x2*z1 + x3*z2 + (x1|group) + (x2|group) + (x3|group), REML = FALSE, data = simdata[[1]], start = list(theta = matrix(c(rep(0, 9)), ncol = 3)))
summary(lme4_fit)
#####################
# stan4bart package #
#####################
stan4bart(y ~ bart(z1*(x1 + x2) + x3*z2 + x4 + x5 + x6 + x7) + (x1|group) + (x2|group) + (x3|group), data = simdata[[1]])
##################
# dbarts package #
##################
dbarts_fit <- rbart_vi(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + z1 + z2,
                       group.by = group,
                       data = simdatasets$simdata_ngroup_30_groupsize_5_icc_0.5_mar_mcar_0_g_0.2[[1]],
                       keepTrees = TRUE,
                       verbose = FALSE)
dbarts_trees <- extract(dbarts_fit, 'trees')
print(dbarts_trees) %>% tail()

dbarts_fit <- bart(x.train = simdata[[1]][,3:11] %>% as.matrix(),
                   y.train = simdata[[1]]$y,
                   ntree = 50,
                   keeptrees = TRUE,
                   nskip = 100,
                   ndpost = 1000)

summary(dbarts_fit)
dbarts_fit$yhat.train.mean
dbarts_fit$ranef.mean

dbarts_treedat <- extractTreeData(dbarts_fit, simdata[[1]])
dbarts_treedat$structure
plotTree(dbarts_treedat, treeNo = 1, iter = 1, plotType = "dendrogram")
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
pred <- tibble(y = simdata[[1]]$y,
               ypred = BART_fit$yhat.train.mean)

ggplot(data = pred, mapping = aes(
  x = y,
  y = ypred
)) +
  geom_point(color = 'royalblue', size = .5) +
  theme_minimal() +
  geom_abline(slope = 1, intercept = 0, color = 'royalblue4') +
  geom_rug(color = 'royalblue', linewidth = .1)

