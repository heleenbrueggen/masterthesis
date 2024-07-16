###################
# Simulation data #
###################
#############
# Libraries #
#############
library(MASS)
library(tidyverse)
library(dplyr)
library(magrittr)
library(psych)
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
################
# Setting path # 
################
path <- "/Volumes/Heleen 480GB/MBART-MICE files/"
###############################
# Simulating multilevel data  #
###############################
###############################################
# Function for calculating the variance of u0 #
###############################################
# Function for calculating the variance for a given correlation
covariance <- function(var1, var2, cor) {
  return(cor * sqrt(var1) * sqrt(var2))
}
# Function for calculating the variance of u0 for specific icc value
var.u0 <- function(t00, phiw, phib, gw, gb, T, sigma2, Sigma, icc) {
  daticc <- icc - ((t(gb) %*% phib %*% gb + t00)/ (t(gw) %*% phiw %*% gw + tr(T %*% Sigma) + t00 + sigma2))

  return(daticc)
}
###########################
# Defining design factors #
###########################
ngroups <- c(30, 50) # Number of groups 
groupsizes <- c(15, 50) # Group sizes 
iccs <- c(.5) # Intraclass correlation coefficient
mar_mcar <- c("mar", "mcar") # Missing data mechanism
miss <- c(0, 50) # Percentage of missing data
g <- c(.5) # Within-group effect size
combinations <- expand.grid(
  ngroup = ngroups,
  groupsize = groupsizes,
  icc = iccs,
  mar_mcar = mar_mcar,
  miss = miss,
  g = g
)
# Saving combinations
write_rds(combinations, file = paste(path, "data/combinations.rds", sep = ""))
###################
# Simulating data #
###################
for (i in seq_len(nrow(combinations))) { # For each combination ...
  # Logging iteration
  cat("Processing iteration:", i, "\n")
  ngroup <- combinations$ngroup[i] # Extracting current number of groups 
  groupsize <- combinations$groupsize[i] # Extracting current group size
  icc <- combinations$icc[i] # Extracting current icc
  # Overall intercept
  g00 <- 10
  # Extracting current within-group effect sizes
  g10 <- combinations$g[i]
  g20 <- combinations$g[i]
  g30 <- combinations$g[i]
  g40 <- combinations$g[i]
  g50 <- combinations$g[i]
  g60 <- combinations$g[i]
  g70 <- combinations$g[i]
  if (icc != 0) { # If icc is not 0 ...
    g01 <- .5
    g02 <- .5
    g11 <- .35
    g21 <- .35
    g32 <- .35
  } else { # If icc is 0 ...
    g01 <- 0
    g02 <- 0
    g11 <- 0
    g21 <- 0
    g32 <- 0
  }

  # Defining matrices for icc calculation (These matrices are explained in more detail in the manuscript)
  phiw <- matrix(c(
    6.25, 2.25, 1.5, 2.55, 1.5, 1.125, 3.3, 0, 0, 0,
    2.25, 9, 1.8, 3.06, 1.8, 1.35, 3.96, 0, 0, 0,
    1.5, 1.8, 4, 2.04, 1.2, .9, 2.64, 0, 0, 0,
    2.25, 3.06, 2.04, 11.56, 2.04, 1.53, 4.488, 0, 0, 0,
    1.5, 1.8, 1.2, 2.04, 4, .9, 2.64, 0, 0, 0,
    1.125, 1.35, 0.9, 1.53, .9, 2.25, 1.98, 0, 0, 0,
    3.3, 3.96, 2.64, 4.488, 2.64, 1.98, 19.36, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 6.25, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 9, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 10.24
  ), nrow = 10, ncol = 10, byrow = TRUE)

  phib <- matrix(c(
    0, 0, 0,
    0, 1, .48,
    0, .48, 2.56
  ), nrow = 3, ncol = 3, byrow = TRUE)

  gw <- c(g10, g20, g30, g40, g50, g60, g70, g11, g21, g32)
  gb <- c(g00, g01, g02)

  T <- matrix(c(
    1, .3, .3, .3, 0, 0, 0, 0,
    .3, 1, .3, .3, 0, 0, 0, 0,
    .3, .3, 1, .3, 0, 0, 0, 0,
    .3, .3, .3, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0
  ), nrow = 8, ncol = 8)

  sigma2 <- 25
  Sigma <- matrix(c(
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 6.25, 2.25, 1.5, 2.55, 1.5, 1.125, 3.3,
    0, 2.25, 9, 1.8, 3.06, 1.8, 1.35, 3.96,
    0, 1.5, 1.8, 4, 2.04, 1.2, .9, 2.64,
    0, 2.25, 3.06, 2.04, 11.56, 2.04, 1.53, 4.488,
    0, 1.5, 1.8, 1.2, 2.04, 4, .9, 2.64,
    0, 1.125, 1.35, 0.9, 1.53, .9, 2.25, 1.98,
    0, 3.3, 3.96, 2.64, 4.488, 2.64, 1.98, 19.36
  ), nrow = 8, ncol = 8)
  
  # Calculating t00, which is the variance of u0
  if (icc != 0) { # If icc is not 0 ...
    t00 <- uniroot(var.u0,
      interval = c(0, 100),
      tol = .00001,
      extendInt = "yes",
      maxiter = 1000,
      phiw = phiw,
      phib = phib,
      gw = gw,
      gb = gb,
      T = T,
      sigma2 = sigma2,
      Sigma = Sigma,
      icc = icc
    )$root
  } else {
    t00 <- 1
  }

  # Simulating data 
  simdata <- replicate(
    n = 1000, # 1000 replications 
    expr = data.frame(
      # individual id
      id = 1:(ngroup * groupsize),
      # group id
      group = rep(1:ngroup, each = groupsize), 
      # residual variance
      eij = rnorm(n = ngroup * groupsize, mean = 0, sd = 5),
      # level-1 predictors
      mvrnorm(
        n = ngroup * groupsize, mu = rep(0, 7),
        Sigma = matrix(c(
          6.25, 2.25, 1.5, 2.55, 1.5, 1.125, 3.3,
          2.25, 9, 1.8, 3.06, 1.8, 1.35, 3.96,
          1.5, 1.8, 4, 2.04, 1.2, .9, 2.64,
          2.25, 3.06, 2.04, 11.56, 2.04, 1.53, 4.488,
          1.5, 1.8, 1.2, 2.04, 4, .9, 2.64,
          1.125, 1.35, 0.9, 1.53, .9, 2.25, 1.98,
          3.3, 3.96, 2.64, 4.488, 2.64, 1.98, 19.36
        ), nrow = 7, ncol = 7, byrow = TRUE)
      ) %>%
        as.data.frame() %>%
        rename(
          x1 = `V1`,
          x2 = `V2`,
          x3 = `V3`,
          x4 = `V4`,
          x5 = `V5`,
          x6 = `V6`,
          x7 = `V7`
        ),
        # level-2 predictors
      mvrnorm(
        n = ngroup, mu = rep(0, 2),
        Sigma = matrix(c(
          1, .48,
          .48, 2.56
        ), nrow = 2, ncol = 2, byrow = TRUE)
      ) %>%
        rep(each = groupsize) %>%
        matrix(ncol = 2) %>%
        as.data.frame() %>%
        rename(
          z1 = `V1`,
          z2 = `V2`
        ),
        # random intercept and slopes
      mvrnorm(
        n = ngroup, mu = rep(0, 4),
        Sigma = matrix(c(
          t00, .3, .3, .3,
          .3, 1, .3, .3,
          .3, .3, 1, .3,
          .3, .3, .3, 1
        ), nrow = 4, ncol = 4, byrow = TRUE)
      ) %>%
        rep(each = groupsize) %>%
        matrix(ncol = 4) %>%
        as.data.frame() %>%
        rename(
          u0 = `V1`,
          u1 = `V2`,
          u2 = `V3`,
          u3 = `V4`,
        )
    ) %>%
      mutate( # Not including random intercept and slopes if icc is 0
        u0 = if (icc != 0) {
          u0
        } else {
          0
        },
        u1 = if (icc != 0) {
          u1
        } else {
          0
        },
        u2 = if (icc != 0) {
          u2
        } else {
          0
        },
        u3 = if (icc != 0) {
          u3
        } else {
          0
        },
        # coefficient generation (including random slopes and cross-level interactions)
        beta0j = g00 + g01 * z1 + g02 * z2 + u0,
        beta1j = g10 + g11 * z1 + u1,
        beta2j = g20 + g21 * z1 + u2,
        beta3j = g30 + g32 * z2 + u3,
        beta4j = g40,
        beta5j = g50,
        beta6j = g60,
        beta7j = g70,
        # generation of dependent variable y
        y = beta0j + beta1j * x1 + beta2j * x2 + beta3j * x3 + beta4j * x4 + beta5j * x5 + beta6j * x6 + beta7j * x7 + eij
      ) %>%
      as_tibble(),
    simplify = FALSE
  )

  # Generate names for each dataset
  name <- paste("simdata",
    colnames(combinations)[1], combinations[i, 1],
    colnames(combinations)[2], combinations[i, 2],
    colnames(combinations)[3], combinations[i, 3],
    colnames(combinations)[4], combinations[i, 4],
    colnames(combinations)[5], combinations[i, 5],
    colnames(combinations)[6], combinations[i, 6],
    sep = "_"
  )
  
  # Saving data in appropriate data folder
  if (combinations$miss[i] != 0) {
    write_rds(simdata, file = paste(path, "data/complete/", name, ".rds", sep = ""))
  } else {
    write_rds(simdata, file = paste(path, "data/nomissing/", name, ".rds", sep = ""))
  }
}
