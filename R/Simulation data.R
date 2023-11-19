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
#####################
# Define parameters #
#####################
ngroup <- 50 # number of clusters
groupsize <- 50 # cluster size
icc <- .5 # intraclass correlation
# fixed effects
g00 <- 1
g01 <- .5
g10 <- 2
g11 <- .3
g20 <- .4
g21 <- 2.7
g30 <- .5
g32 <- 1.8
g40 <- .2
g50 <- 3
g60 <- 1.3
g70 <- .4
simdata <- simdata(nsim = 1000,
                   ngroup = 50,
                   groupsize = 50,
                   icc = .5,
                   g00 = 1, g01 = .5, g10 = 2, g11 = .3, g20 = .4, g21 = 2.7, g30 = .5, g32 = 1.8, g40 = .2, g50 = 3, g60 = 1.3, g70 = .4,
                   vare = 1,
                   varx1 = 2, varx2 = 3.1, varx3 = 1.5, varx4 = 4, varx5 = 6, varx6 = 2.5, varx7 = 2,
                   varz1 = 1.5, varz2 = 5,
                   varu1 = 1, varu2 = 1, varu3 = 1, varu4 = 1, varu5 = 1, varu6 = 1,
                   meane = 0, meanx1 = 0, meanx2 = 0, meanx3 = 0, meanx4 = 0, meanx5 = 0, meanx6 = 0, meanx7 = 0,
                   meanu0 = 0, meanu1 = 0, meanu2 = 0, meanu3 = 0, meanu4 = 0, meanu5 = 0, meanu6 = 0)
#############################
# Generate simulation data  #
#############################
#####
# Simulation data with only one x and z variable
#####
simdata <- lapply(1:1000, function(i) {
  data <- tibble(
    id = 1:(ngroup * groupsize),
    group = rep(1:ngroup, each = groupsize),
    eij = rnorm(n = ngroup * groupsize, mean = 0, sd = 5),
    x1 = rnorm(n = ngroup * groupsize, mean = 0, sd = 2.5),
    z = rep(rnorm(n = ngroup, mean = 0, sd = 1), each = groupsize),
    u1 = rep(rnorm(n = ngroup, mean = 0, sd = 1), each = groupsize),
    u0 = rep(rnorm(n = ngroup,
                   mean = 0,
                   sd = uniroot(function(varu0, g00, g01, z, g10, g11, u1, x1, eij, icc) {

      daticc <- (icc - ((var(g00 + g01 * z) + varu0) / (var(g00 + g01 * z) + varu0 + var((g10 + g11 * z + u1)*x1) + var(eij))))

      return(daticc)
    }, interval = c(0, 100),
    tol = .0001,
    extendInt = 'yes',
    maxiter = 1000, g00 = g00, g01 = g01, z = z, g10 = g10, g11 = g11, u1 = u1, x1 = x1, eij = eij, icc = icc)$root %>% as.numeric() %>% sqrt()), each = groupsize),
    beta0j = g00 + g01 * z + u0,
    beta1j = g10 + g11 * z + u1,
    y = beta0j + beta1j * x1 + eij)
}
)
#####
# Simulation like in thesis proposal
#####
simdata <- lapply(1:1000, function(i) {
  data <- tibble(
    id = 1:(ngroup * groupsize),
    group = rep(1:ngroup, each = groupsize),
    eij = rnorm(n = ngroup * groupsize, mean = 0, sd = 5),
    x1 = rnorm(n = ngroup * groupsize, mean = 0, sd = 2.5),
    x2 = rnorm(n = ngroup * groupsize, mean = 0, sd = 3),
    z = rep(rnorm(n = ngroup, mean = 0, sd = 1), each = groupsize),
    u2 = rep(rnorm(n = ngroup, mean = 0, sd = 1), each = groupsize),
    u1 = rep(rnorm(n = ngroup, mean = 0, sd = 1), each = groupsize),
    u0 = rep(rnorm(n = ngroup,
                   mean = 0,
                   sd = uniroot(function(varu0, g00, g01, z, g10, g11, u1, x1, g20, u2, eij, icc) {

                     daticc <- (icc - ((var(g00 + g01 * z) + varu0) / (var(g00 + g01 * z) + varu0 + var((g10 + g11 * z + u1)*x1) + var((g20 + u2)*x2) + var(eij))))

                     return(daticc)
                   }, interval = c(0, 100),
                   tol = .0001,
                   extendInt = 'yes',
                   maxiter = 1000, g00 = g00, g01 = g01, z = z, g10 = g10, g11 = g11, u1 = u1, x1 = x1, g20 = g20, u2 = u2, eij = eij, icc = icc)$root %>% as.numeric() %>% sqrt()), each = groupsize),
    beta0j = g00 + g01 * z + u0,
    beta1j = g10 + g11 * z + u1,
    y = beta0j + beta1j * x1 + eij)
}
)

simdata <- replicate(n = 1000,
                     expr = tibble(
                       id = 1:(ngroup * groupsize),
                       group = rep(1:ngroup, each = groupsize),
                       eij = rnorm(n = ngroup * groupsize, mean = 0, sd = 5),
                       x1 = rnorm(n = ngroup * groupsize, mean = 0, sd = 2.5),
                       x2 = rnorm(n = ngroup * groupsize, mean = 0, sd = 3),
                       z = rep(rnorm(n = ngroup, mean = 0, sd = 1), each = groupsize),
                       u2 = rep(rnorm(n = ngroup, mean = 0, sd = 1), each = groupsize),
                       u1 = rep(rnorm(n = ngroup, mean = 0, sd = 1), each = groupsize),
                       u0 = rep(rnorm(n = ngroup,
                                      mean = 0,
                                      sd = uniroot(function(varu0, g00, g01, z, g10, g11, u1, x1, g20, u2, eij, icc) {

                                        daticc <- (icc - ((var(g00 + g01 * z) + varu0) / (var(g00 + g01 * z) + varu0 + var((g10 + g11 * z + u1)*x1) + var((g20 + u2)*x2) + var(eij))))

                                        return(daticc)
                                      }, interval = c(0, 100),
                                      tol = .0001,
                                      extendInt = 'yes',
                                      maxiter = 1000, g00 = g00, g01 = g01, z = z, g10 = g10, g11 = g11, u1 = u1, x1 = x1, g20 = g20, u2 = u2, eij = eij, icc = icc)$root %>% as.numeric() %>% sqrt()), each = groupsize),
                       beta0j = g00 + g01 * z + u0,
                       beta1j = g10 + g11 * z + u1,
                       beta2j = g20 + u2,
                       y = beta0j + beta1j * x1 + beta2j * x2 + eij),
                     simplify = FALSE)
#####
# Simulation complicated model with 7 variables first level, 3 second level
#####
simdata <- function (nsim = 1000,
                     ngroup,
                     groupsize,
                     icc,
                     ...) {
  # Input values:
  # nsim = number of desired simulated data sets
  # ngroup = number of desired groups
  # groupsize = desired size of the groups
  # icc = desired icc in the data set
  # ... = additional arguments pertaining to the fixed effects, means and variances

  # Extracting additional arguments using list(...)
  args <- list(...)

  # Simulated data generation
  data <- replicate(n = nsim,
                     expr = tibble(
                       # individual id
                       id = 1:(ngroup * groupsize),
                       # group id
                       group = rep(1:ngroup, each = groupsize),
                       # residual variance
                       eij = rnorm(n = ngroup * groupsize, mean = args$meane, sd = sqrt(args$vare)),
                       # level 1 variables
                       x1 = rnorm(n = ngroup * groupsize, mean = args$meanx1, sd = sqrt(args$varx1)),
                       x2 = rnorm(n = ngroup * groupsize, mean = args$meanx2, sd = sqrt(args$varx2)),
                       x3 = rnorm(n = ngroup * groupsize, mean = args$meanx3, sd = sqrt(args$varx3)),
                       x4 = rnorm(n = ngroup * groupsize, mean = args$meanx4, sd = sqrt(args$varx4)),
                       x5 = rnorm(n = ngroup * groupsize, mean = args$meanx5, sd = sqrt(args$varx5)),
                       x6 = rnorm(n = ngroup * groupsize, mean = args$meanx6, sd = sqrt(args$varx6)),
                       x7 = rnorm(n = ngroup * groupsize, mean = args$meanx7, sd = sqrt(args$varx7)),
                       # level 2 variables
                       z1 = rep(rnorm(n = ngroup, mean = args$meanz1, sd = sqrt(args$varz1)), each = groupsize),
                       z2 = rep(rnorm(n = ngroup, mean = args$meanz2, sd = sqrt(args$varz2)), each = groupsize),
                       # random slopes
                       u6 = rep(rnorm(n = ngroup, mean = args$meanu6, sd = sqrt(args$varu6)), each = groupsize),
                       u5 = rep(rnorm(n = ngroup, mean = args$meanu5, sd = sqrt(args$varu5)), each = groupsize),
                       u4 = rep(rnorm(n = ngroup, mean = args$meanu4, sd = sqrt(args$varu4)), each = groupsize),
                       u3 = rep(rnorm(n = ngroup, mean = args$meanu3, sd = sqrt(args$varu3)), each = groupsize),
                       u2 = rep(rnorm(n = ngroup, mean = args$meanu2, sd = sqrt(args$varu2)), each = groupsize),
                       u1 = rep(rnorm(n = ngroup, mean = args$meanu1, sd = sqrt(args$varu1)), each = groupsize),
                       u0 = rep(rnorm(n = ngroup,
                                      mean = args$meanu0,
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
                                      maxiter = 1000, g00=args$g00, g01=args$g01, g11=args$g11, g20=args$g20, g21=args$g21, g30=args$g30, g32=args$g32, g40=args$g40, g50=args$g50, g60=args$g60, g70=args$g70, z1=z1, z2=z2, x1=x1, x2=x2, x3=x3, x4=x4, x5=x5, x6=x6, x7=x7, eij=eij, icc=args$icc)$root %>% as.numeric() %>% sqrt()), each = groupsize),
                       # coefficient generation (including random slopes and cross-level interactions)
                       beta0j = args$g00 + args$g01 * z1 + u0,
                       beta1j = args$g10 + args$g11 * z1 + u1,
                       beta2j = args$g20 + args$g21 * z1 + u2,
                       beta3j = args$g30 + args$g32 * z2 + u3,
                       beta4j = args$g40 + u4,
                       beta5j = args$g50 + u5,
                       beta6j = args$g60 + u6,
                       beta7j = args$g70,
                       # generation of dependent variable y
                       y = beta0j + beta1j * x1 + beta2j * x2 + beta3j * x3 + beta4j * x4 + beta5j * x5 + beta6j * x6 * beta7j * x7 + eij) %>%
                       # taking out terms that are only used for model generation
                       select(-u0, -u1, -u2, -u3, -u4, -u5, -u6, -eij, -beta0j, -beta1j, -beta2j, -beta3j, -beta4j, -beta5j, -beta6j, -beta7j),
                     simplify = FALSE)

  # Output
  return(data) # list of simulated data sets
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

