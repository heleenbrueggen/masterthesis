#Stan4bart

# Installing package
# devtools::install_github("vdorie/stan4bart")

# Libraries
library(stan4bart)
library(BART)
library(lme4)
library(dbarts)

################
# lme4 package #
################
lme4_fit <- lmer(y ~ 1 + x1 + x2 + x3 + x4 + x5 + x6 + x7 + z1 + z2 + x1*z1 + x2*z1 + x3*z2 + (x1|group) + (x2|group) + (x3|group), REML = FALSE, data = simdata[[1]])
summary(lme4_fit)
#####################
# stan4bart package #
#####################
stan4bart(y ~ bart(z1*(x1 + x2) + x3*z2 + x4 + x5 + x6 + x7) + (x1|group) + (x2|group) + (x3|group), data = simdata[[1]])
##################
# dbarts package #
##################
dbarts_fit <- rbart_vi(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + z1 + z2, group.by = group, data = simdata[[1]])
summary(dbarts_fit)
dbarts_fit$yhat.train.mean
dbarts_fit$ranef.mean
################
# BART package #
################
BART_fit <- wbart(x.train = simdata[[1]][,3:11] %>% as.matrix(), y.train = simdata[[1]]$y)
BART_fit$yhat.train.mean # train data fits
