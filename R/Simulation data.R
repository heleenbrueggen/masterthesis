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
library(readr)
################
# Setting seed #
################
set.seed(123)
###############################
# Simulating multilevel data  #
###############################
###############################################
# Function for calculating the variance of u0 #
###############################################
var.u0 <- function(varu0,
                   g00,
                   g01, g02, g11, g21, g32,
                   g10, g20, g30, g40, g50, g60, g70,
                   z1, z2,
                   x1, x2, x3, x4, x5, x6, x7,
                   u1, u2, u3, u4, u5, u6,
                   eij,
                   icc) {
  daticc <- icc - ((var(g00 + g01 * z1 + g02 * z2) + varu0) / (var(g00 + g01 * z1 + g02 * z2) + varu0 +
    var((g10 + g11 * z1 + u1) * x1) +
    var((g20 + g21 * z1 + u2) * x2) +
    var((g30 + g32 * z2 + u3) * x3) +
    var((g40 + u4) * x4) +
    var((g50 + u5) * x5) +
    var((g60 + u6) * x6) +
    var(g70 * x7) + var(eij)))
  return(daticc)
}
#######################
# Defining parameters #
#######################
ngroups <- c(30, 50)
groupsizes <- c(5, 15, 35, 50)
iccs <- c(0, .05, .3, .5)
mar_mcar <- c("mar", "mcar")
miss <- c(0, 25, 50)
g <- c(.2, .5, .8)
combinations <- expand.grid(
  ngroup = ngroups,
  groupsize = groupsizes,
  icc = iccs,
  mar_mcar = mar_mcar,
  miss = miss,
  g = g
)
###################
# Simulating data #
###################
for (i in seq_len(nrow(combinations))) {
  # Logging iteration
  cat("Processing iteration:", i, "\n")
  ngroup <- combinations$ngroup[i]
  groupsize <- combinations$groupsize[i]
  icc <- combinations$icc[i]
  # Overall intercept
  g00 <- 10
  # Individual effects
  g10 <- combinations$g[i]
  g20 <- combinations$g[i]
  g30 <- combinations$g[i]
  g40 <- combinations$g[i]
  g50 <- combinations$g[i]
  g60 <- combinations$g[i]
  g70 <- combinations$g[i]
  if (icc != 0) {
    g01 <- .5
    g02 <- .5
    g11 <- .35
    g21 <- .35
    g32 <- .35
  } else {
    g01 <- 0
    g02 <- 0
    g11 <- 0
    g21 <- 0
    g32 <- 0
  }

  simdata <- replicate(
    n = 10,
    expr = tibble(
      # individual id
      id = 1:(ngroup * groupsize),
      # group id
      group = rep(1:ngroup, each = groupsize),
      # residual variance
      eij = rnorm(n = ngroup * groupsize, mean = 0, sd = 5),
      # level 1 variables
      x1 = rnorm(n = ngroup * groupsize, mean = 0, sd = 2.5),
      x2 = rnorm(n = ngroup * groupsize, mean = 0, sd = 3),
      x3 = rnorm(n = ngroup * groupsize, mean = 0, sd = 2),
      x4 = rnorm(n = ngroup * groupsize, mean = 0, sd = 3.4),
      x5 = rnorm(n = ngroup * groupsize, mean = 0, sd = 2),
      x6 = rnorm(n = ngroup * groupsize, mean = 0, sd = 1.5),
      x7 = rnorm(n = ngroup * groupsize, mean = 0, sd = 4.4),
      # level 2 variables
      z1 = rep(rnorm(n = ngroup, mean = 0, sd = 1), each = groupsize),
      z2 = rep(rnorm(n = ngroup, mean = 0, sd = 1.6), each = groupsize),
      # random slopes
      u6 = if (icc != 0) {
        rep(rnorm(n = ngroup, mean = 0, sd = 1), each = groupsize)
      } else {
        rep(0, (ngroup * groupsize))
      },
      u5 = if (icc != 0) {
        rep(rnorm(n = ngroup, mean = 0, sd = 1), each = groupsize)
      } else {
        rep(0, (ngroup * groupsize))
      },
      u4 = if (icc != 0) {
        rep(rnorm(n = ngroup, mean = 0, sd = 1), each = groupsize)
      } else {
        rep(0, (ngroup * groupsize))
      },
      u3 = if (icc != 0) {
        rep(rnorm(n = ngroup, mean = 0, sd = 1), each = groupsize)
      } else {
        rep(0, (ngroup * groupsize))
      },
      u2 = if (icc != 0) {
        rep(rnorm(n = ngroup, mean = 0, sd = 1), each = groupsize)
      } else {
        rep(0, (ngroup * groupsize))
      },
      u1 = if (icc != 0) {
        rep(rnorm(n = ngroup, mean = 0, sd = 1), each = groupsize)
      } else {
        rep(0, (ngroup * groupsize))
      },
      u0 = if (icc != 0) {
        rep(rnorm(
          n = ngroup,
          mean = 0,
          sd = uniroot(var.u0,
            interval = c(0, 100),
            tol = .000001,
            extendInt = "yes",
            maxiter = 1000,
            g00 = g00,
            g01 = g01, g02 = g02, g11 = g11, g21 = g21, g32 = g32,
            g10 = g10, g20 = g20, g30 = g30, g40 = g40, g50 = g50, g60 = g60, g70 = g70,
            z1 = z1, z2 = z2,
            x1 = x1, x2 = x2, x3 = x3, x4 = x4, x5 = x5, x6 = x6, x7 = x7,
            u1 = u1, u2 = u2, u3 = u3, u4 = u4, u5 = u5, u6 = u6,
            eij = eij,
            icc = icc
          )$root %>% sqrt()
        ), each = groupsize)
      } else {
        rep(0, (ngroup * groupsize))
      },
      # coefficient generation (including random slopes and cross-level interactions)
      beta0j = g00 + g01 * z1 + g02 * z2 + u0,
      beta1j = g10 + g11 * z1 + u1,
      beta2j = g20 + g21 * z1 + u2,
      beta3j = g30 + g32 * z2 + u3,
      beta4j = g40 + u4,
      beta5j = g50 + u5,
      beta6j = g60 + u6,
      beta7j = g70,
      # generation of dependent variable y
      y = beta0j + beta1j * x1 + beta2j * x2 + beta3j * x3 + beta4j * x4 + beta5j * x5 + beta6j * x6 * beta7j * x7 + eij
    ) %>%
      # taking out terms that are only used for model generation
      select(-u0, -u1, -u2, -u3, -u4, -u5, -u6, -eij, -beta0j, -beta1j, -beta2j, -beta3j, -beta4j, -beta5j, -beta6j, -beta7j),
    simplify = FALSE
  )
  
  name <- paste("simdata",
    colnames(combinations)[1], combinations[i, 1],
    colnames(combinations)[2], combinations[i, 2],
    colnames(combinations)[3], combinations[i, 3],
    colnames(combinations)[4], combinations[i, 4],
    colnames(combinations)[5], combinations[i, 5],
    colnames(combinations)[6], combinations[i, 6],
    sep = "_"
  )
  assign(name, simdata)
  # Saving data in data folder
  write_rds(simdata, file = paste("data/", name, ".rds", sep = ""))
}
#############
# Load data #
#############
simdatasets <- list()
for (i in seq_len(nrow(combinations))) {
  name <- paste("simdata",
    colnames(combinations)[1], combinations[i, 1],
    colnames(combinations)[2], combinations[i, 2],
    colnames(combinations)[3], combinations[i, 3],
    colnames(combinations)[4], combinations[i, 4],
    colnames(combinations)[5], combinations[i, 5],
    colnames(combinations)[6], combinations[i, 6],
    sep = "_"
  )
  simdatasets[[i]] <- read_rds(paste("data/", name, ".rds", sep = ""))
}
#############################
# Storing names of datasets #
#############################
names <- rep(NA, 576)
for (i in seq_len(nrow(combinations))) {
  names[i] <- paste("simdata",
    colnames(combinations)[1], combinations[i, 1],
    colnames(combinations)[2], combinations[i, 2],
    colnames(combinations)[3], combinations[i, 3],
    colnames(combinations)[4], combinations[i, 4],
    colnames(combinations)[5], combinations[i, 5],
    colnames(combinations)[6], combinations[i, 6],
    sep = "_"
  )
}
####################
# Check simulation #
####################
#########################################
# Check for population generation model #
#########################################
# Checking population generation model over all simulated data sets
simdatasets %>%
  map(~.x %$%
        lmer(y ~ 1 + x1 + x2 + x3 + x4 + x5 + x6 + x7 + z1 + z2 + x1 * z1 + x2 * z1 + x3 * z2 + (x1 | group) + (x2 | group) + (x3 | group), REML = FALSE) %>%
        summary %>%
        coefficients %>%
        .[, 1]) %>%
  Reduce("+", .) / length(simdatasets)
################
# Checking ICC #
################
# Checking ICC values of all datasets
iccvalues <- rep(0, 576)
for (i in 1:288) {
  iccvalues[i] <- simdatasets[[i]] %>%
    map(~ .x %$%
      lmer(y ~ 1 + (1 | group), REML = FALSE) %>%
      summ() %>%
      .$gvars %>%
      .[3]) %>%
    as.numeric() %>%
    Reduce("+", .) / length(simdatasets[[i]])
}
# Check if they are the same as specified
cbind(iccvalues %>% round(digits = 3), combinations)
############
# Appendix #
############
