###################
# Simulation data #
###################
#############
# Libraries #
#############
library(tidyverse)
library(dplyr)
library(magrittr)
library(psych)
library(mvtnorm)
library(lme4)
library(lmerTest)
library(jtools)
library(purrr)
library(tibble)
################
# Setting seed #
################
set.seed(123)
#############################
# Generate simulation data  #
#############################
# Simulation function complicated model with 7 variables first level, 3 second level
simdata <- function (nsim = 1000,
                     ngroup,
                     groupsize,
                     icc,
                     g00 = 1, g01 = 1, g10 = 1, g11 = 1, g20 = 1, g21 = 1, g30 = 1, g32 = 1, g40 = 1, g50 = 1, g60 = 1, g70 = 1,
                     vare = 1, varx1 = 1, varx2 = 1, varx3 = 1, varx4 = 1, varx5 = 1, varx6 = 1, varx7 = 1, varz1 = 1, varz2 = 1, varu1 = 1, varu2 = 1, varu3 = 1, varu4 = 1, varu5 = 1, varu6 = 1,
                     meane = 0, meanx1 = 0, meanx2 = 0, meanx3 = 0, meanx4 = 0, meanx5 = 0, meanx6 = 0, meanx7 = 0, meanz1 = 0, meanz2 = 0, meanu0 = 0, meanu1 = 0, meanu2 = 0, meanu3 = 0, meanu4 = 0, meanu5 = 0, meanu6 = 0) {
  # Input values:
  # nsim = number of desired simulated data sets
  # ngroup = number of desired groups
  # groupsize = desired size of the groups
  # icc = desired icc in the data set
  # g = fixed effects,
  # mean = means,
  # var = variances

  # Simulated data generation
  data <- replicate(n = nsim,
                     expr = tibble(
                       # individual id
                       id = 1:(ngroup * groupsize),
                       # group id
                       group = rep(1:ngroup, each = groupsize),
                       # residual variance
                       eij = rnorm(n = ngroup * groupsize, mean = meane, sd = sqrt(vare)),
                       # level 1 variables
                       x1 = rnorm(n = ngroup * groupsize, mean = meanx1, sd = sqrt(varx1)),
                       x2 = rnorm(n = ngroup * groupsize, mean = meanx2, sd = sqrt(varx2)),
                       x3 = rnorm(n = ngroup * groupsize, mean = meanx3, sd = sqrt(varx3)),
                       x4 = rnorm(n = ngroup * groupsize, mean = meanx4, sd = sqrt(varx4)),
                       x5 = rnorm(n = ngroup * groupsize, mean = meanx5, sd = sqrt(varx5)),
                       x6 = rnorm(n = ngroup * groupsize, mean = meanx6, sd = sqrt(varx6)),
                       x7 = rnorm(n = ngroup * groupsize, mean = meanx7, sd = sqrt(varx7)),
                       # level 2 variables
                       z1 = rep(rnorm(n = ngroup, mean = meanz1, sd = sqrt(varz1)), each = groupsize),
                       z2 = rep(rnorm(n = ngroup, mean = meanz2, sd = sqrt(varz2)), each = groupsize),
                       # random slopes
                       u6 = rep(rnorm(n = ngroup, mean = meanu6, sd = sqrt(varu6)), each = groupsize),
                       u5 = rep(rnorm(n = ngroup, mean = meanu5, sd = sqrt(varu5)), each = groupsize),
                       u4 = rep(rnorm(n = ngroup, mean = meanu4, sd = sqrt(varu4)), each = groupsize),
                       u3 = rep(rnorm(n = ngroup, mean = meanu3, sd = sqrt(varu3)), each = groupsize),
                       u2 = rep(rnorm(n = ngroup, mean = meanu2, sd = sqrt(varu2)), each = groupsize),
                       u1 = rep(rnorm(n = ngroup, mean = meanu1, sd = sqrt(varu1)), each = groupsize),
                       u0 = rep(rnorm(n = ngroup,
                                      mean = meanu0,
                                      sd = uniroot(function(varu0, g00, g01, g11, g20, g21, g30, g32, g40, g50, g60, g70, z1, z2, x1, x2, x3, x4, x5, x6, x7, eij, icc) {

                                        daticc <- (icc - ((var(g00 + g01 * z1) + varu0) / (var(g00 + g01 * z1) + varu0 +
                                                                                             var((g10 + g11 * z1 + u1)*x1) +
                                                                                             var((g20 + g21 * z1 + u2)*x2) +
                                                                                             var((g30 + g32 * z2 + u3)*x3) +
                                                                                             var((g40 + u4)*x4) +
                                                                                             var((g50 + u5)*x5) +
                                                                                             var((g60 + u6)*x6) +
                                                                                             var(g70*x7) + var(eij))))

                                        return(daticc)
                                      }, interval = c(0, 100),
                                      tol = .0001,
                                      extendInt = 'yes',
                                      maxiter = 1000, g00=g00, g01=g01, g11=g11, g20=g20, g21=g21, g30=g30, g32=g32, g40=g40, g50=g50, g60=g60, g70=g70, z1=z1, z2=z2, x1=x1, x2=x2, x3=x3, x4=x4, x5=x5, x6=x6, x7=x7, eij=eij, icc=icc)$root %>% as.numeric() %>% sqrt()), each = groupsize),
                       # coefficient generation (including random slopes and cross-level interactions)
                       beta0j = g00 + g01 * z1 + u0,
                       beta1j = g10 + g11 * z1 + u1,
                       beta2j = g20 + g21 * z1 + u2,
                       beta3j = g30 + g32 * z2 + u3,
                       beta4j = g40 + u4,
                       beta5j = g50 + u5,
                       beta6j = g60 + u6,
                       beta7j = g70,
                       # generation of dependent variable y
                       y = beta0j + beta1j * x1 + beta2j * x2 + beta3j * x3 + beta4j * x4 + beta5j * x5 + beta6j * x6 * beta7j * x7 + eij) %>%
                       # taking out terms that are only used for model generation
                       select(-u0, -u1, -u2, -u3, -u4, -u5, -u6, -eij, -beta0j, -beta1j, -beta2j, -beta3j, -beta4j, -beta5j, -beta6j, -beta7j),
                     simplify = FALSE)

  # Output
  return(data) # list of simulated data sets
}
############################
# Simulating all data sets #
############################
ngroups <- c(30, 50)
groupsizes <- c(5, 15, 35, 50)
iccs <- c(0, .05, .3, .5)
mar_mcar <- c(0, 25, 50)
g <- c(.2, .5, .8)
combinations <- expand.grid(ngroup = ngroups,
                            groupsize = groupsizes,
                            icc = iccs,
                            mar_mcar = mar_mcar,
                            g = g)

simdatasets <- list()
for (i in 1:nrow(combinations)) {
  ngroup <- combinations$ngroup[i]
  groupsize <- combinations$groupsize[i]
  icc <- combinations$icc[i]
  g <- combinations$g[i]

  simdata <- simdata(nsim = 1000, ngroup = ngroup, groupsize = groupsize, icc = icc,
                     g10 = g, g20 = g, g30 = g, g40 = g, g50 = g, g60 = g, g70 = g)

  name <- paste("simdata_", i, sep = "")
  simdatasets[[name]] <- simdata
}
####################
# Check simulation #
####################
#########################################
# Check for population generation model #
#########################################
# Checking population generation model over all simulated data sets
simdata %>%
  map(~.x %$%
        lmer(y ~ 1 + x1 + x2 + x3 + x4 + x5 + x6 + x7 + z1 + z2 + x1*z1 + x2*z1 + x3*z2 + (x1|group) + (x2|group) + (x3|group), REML = FALSE) %>%
        summary %>%
        coefficients %>%
        .[,1]) %>%
  Reduce("+", .) / length(simdata)

################
# Checking ICC #
################
# Generating function for ICC calculation
iccfunction <- function (data) {
  icc <- var(data$beta0j)/var(data$y)
  return(icc)
}
# Checking ICC over all simulated data sets
iccvalue <- rep(0, 1000)
for(i in 1:1000) {
  iccvalue[i] <- iccfunction(simdata[[i]])
}
mean(iccvalue)

simdata %>%
  map(~.x %$%
        lmer(y ~ 1 + (1|group), REML = FALSE) %>%
        summ %>%
        .$gvars %>%
        .[3]) %>% as.numeric() %>%
  Reduce("+", .) / length(simdata)
############
# Appendix #
############
simdata %>%
  map(~.x %$%
        var(x1)) %>%
  Reduce("+", .) / length(simdata)

