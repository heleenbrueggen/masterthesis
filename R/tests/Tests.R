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
#######################
# Defining parameters #
#######################
ngroups <- c(100)
groupsizes <- c(100)
iccs <- c(.5)
mar_mcar <- c("mcar")
miss <- c(50)
g <- c(.5)
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

  # Defining matrices for icc calculation
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

  if (icc != 0) {
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

  simdata <- replicate(
    n = 1000,
    expr = data.frame(
      # individual id
      id = 1:(ngroup * groupsize),
      # group id
      group = rep(1:ngroup, each = groupsize),
      # residual variance
      eij = rnorm(n = ngroup * groupsize, mean = 0, sd = 5),
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
      mutate( 
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
      # taking out terms that are only used for model generation
      # select(-u0, -u1, -u2, -u3, -eij, -beta0j, -beta1j, -beta2j, -beta3j, -beta4j, -beta5j, -beta6j, -beta7j) %>%
      as_tibble(),
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
  
  # Saving data in appropriate data folder
#   if (combinations$miss[i] != 0) {
#     write_rds(simdata, file = paste("data/complete/", name, ".rds", sep = ""))
#   } else {
#     write_rds(simdata, file = paste("data/nomissing/", name, ".rds", sep = ""))
#   }
}
write_rds(simdata, file = "data/test/simdata_test.rds")
#############################
# Storing names of datasets #
#############################
names <- rep(NA, 576)
for (i in seq_len(nrow(combinations))) {
  names[i] <- paste(
    colnames(combinations)[1], combinations[i, 1],
    colnames(combinations)[2], combinations[i, 2],
    colnames(combinations)[3], combinations[i, 3],
    colnames(combinations)[4], combinations[i, 4],
    colnames(combinations)[5], combinations[i, 5],
    colnames(combinations)[6], combinations[i, 6],
    sep = "_"
  )
}
#############
# Load data #
#############
simdatasets <- list()
for (i in seq_len(nrow(combinations))) {
  if (combinations$miss[i] != 0) {
    simdatasets[[i]] <- read_rds(paste("data/complete/simdata_", names[i], ".rds", sep = ""))
  } else {
    simdatasets[[i]] <- read_rds(paste("data/nomissing/simdata_", names[i], ".rds", sep = ""))
  }
}
####################
# Check simulation #
####################
library(lme4)
library(lmerTest)
library(magrittr)
simdata <- read_rds("data/test/simdata_test.rds")

simdata2 <- do.call(rbind, simdata)

psych::describe(simdata[[2]])

test <- lmer(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + z1 + z2 + x1 * z1 + x2 * z1 + x3 * z2 + (1 + x1 + x2 + x3 | group),
    REML = TRUE,
    control = lmerControl(optimizer = "bobyqa"), data = simdata[[2]]
)
summary(test)$coefficients[,1]
broom.mixed::tidy(test, effects = c("fixed", "ran_pars"), scales = c(NA, "vcov"))
jtools::summ(test)

lm(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + z1 + z2 + x1 * z1 + x2 * z1 + x3 * z2, data = simdata[[2]]) %>%
    broom.mixed::tidy()
test
performance::check_model(test)

test2 <- nlme::lme(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + z1 + z2 + x1 * z1 + x2 * z1 + x3 * z2, random = ~1 + x1 + x2 + x3 | group, data = simdata[[2]])

performance::check_model(test2)

test3 <- map(simdata, ~ .x %$% lmer(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + z1 + z2 + x1 * z1 + x2 * z1 + x3 * z2 + (1 + x1 + x2 + x3 | group),
    REML = TRUE,
    control = lmerControl(optimizer = "bobyqa"), data = .
) %>% summary() %$% coefficients[,1]) # %>%
    # Reduce("+", .) / length(simdata)
write_rds(test3, file = "data/test/test3.rds")
test3 <- read_rds("data/test/test3.rds")

str(simdata[[2]])
test4 <- simdata[[2]] %>% mutate(
    g01z1g02z2 = beta0j - u0 - 10,
    z1z2 = g01z1g02z2 / (z1 + z2),
    g01z1 = .5 * z1,
    g02z2 = .5 * z2,
    beta0 = 10 + .5 * z1 + .5 * z2 + u0
)
test4[, c("beta0j", "u0", "g01z1g02z2", "z1z2", "g01z1", "g02z2", "beta0")] %>% View()


simdatasets_nomiss <- read_rds(paste("/Volumes/Heleen\ 480GB/Master\ thesis/data/nomissing/simdata_ngroup_50_groupsize_50_icc_0.5_mar_mcar_mar_miss_0_g_0.5.rds", sep = ""))[1:100]
test5 <- lmer(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + z1 + z2 + x1 * z1 + x2 * z1 + x3 * z2 + (1 + x1 + x2 + x3 | group),
    REML = TRUE,
    control = lmerControl(optimizer = "bobyqa"),
    data = simdatasets_nomiss[[1]]
)
summary(test5)
jtools::summ(test5)

simdatasets_nomiss2 <- read_rds(paste("/Volumes/Heleen\ 480GB/Master\ thesis/data/nomissing/simdata_ngroup_30_groupsize_15_icc_0.5_mar_mcar_mar_miss_0_g_0.5.rds", sep = ""))[1:100]
test6 <- lmer(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + z1 + z2 + x1 * z1 + x2 * z1 + x3 * z2 + (1 + x1 + x2 + x3 | group),
    REML = TRUE,
    control = lmerControl(optimizer = "bobyqa"),
    data = simdatasets_nomiss2[[1]]
)
jtools::summ(test6)

summary(test6)
test7 <- map(simdatasets_nomiss2, ~ .x %$% lmer(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + z1 + z2 + x1 * z1 + x2 * z1 + x3 * z2 + (1 + x1 + x2 + x3 | group),
    REML = TRUE,
    control = lmerControl(optimizer = "bobyqa"), data = .
) %>% summary() %$% coefficients[,1]) %>%
    Reduce("+", .) / length(simdatasets_nomiss2)

library(lavaan)
model <- "
    level: 1
        y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x1:z1 + x2:z1 + x3:z2

    level: 2 
        y ~ z1 + z2
"
fit <- sem(model, data = simdatasets_nomiss2[[1]], cluster = "group", optim.method = "em")
summary(fit)
################
# Checking ICC #
################
# Checking ICC values of all datasets
iccvalues <- rep(0, 576)
for (i in 1:576) {
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
cbind(iccvalues %>% round(digits = 3), combinations$icc)



var.u0 <- function(t00, phiw, phib, gw, gb, T, sigma2, Sigma, icc) {
    daticc <- icc - ((t(gb) %*% phib %*% gb + t00) / (t(gw) %*% phiw %*% gw + tr(T %*% Sigma) + t00 + sigma2))

    return(daticc)
}
# Defining matrices for icc calculation
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

gw <- c(.5, .5, .5, .5, .5, .5, .5, .35, .35, .35)
gb <- c(10, .5, .5)

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
    icc = .5
)$root
