# no scientific notation
options(scipen = 999)

# remotes::install_github('blimp-stats/rblimp')
library(fdir)
library(mlmpower)
library(mvtnorm)
library(rblimp)
library(lmerTest)
library(varTestnlme)
library(ggplot2)
# library(RLRsim)

#################################################################
# POPULATION MODEL:
#
# Yij = (B0 + b1j) + (B1 + b1j)(X1ij.w) + (B2 + b1j)(X2ij.w) + (B3)(X1j.b) + (B4)(X2j.b) + eij
#
# The R code invokes the following constraints: B1 = B3 and B2 = B4
#################################################################

# use mlmpower package to solve for MLM parameters given Rights & Sterba effect sizes (see Enders, Keller, & Woller in Psych Methods)
ranslope <- (
  effect_size(
    icc          = c(.50),
    within       = .135,
    random_slope = .05
  )
  + outcome("y", mean = 50, sd = 10)
    + within_predictor("x1", weight = .5)
    + within_predictor("x2", weight = .5)
    + random_slope("x1", weight = .5)
    + random_slope("x2", weight = .5)
    + correlations(randeff = 0)
)

params_ranslope <- summary(ranslope)

#################################################################
# function to generate 1 data set
#################################################################

gen1data <- function(params, numL2, numperL2) {
  # params = params_null
  # numL2 = 20
  # numperL2 = 10

  ################## extract parameters ##################

  # focal model parameters
  betas <- params$gammas # coefficients
  cov_b0b1 <- params$tau # covariance matrix of random intercepts and slopes
  var_e <- params$var_e # within-cluster residual variance

  # predictor parameters
  mean_X <- params$mean_X # within-cluster predictor mean
  mean_Z <- params$mean_Z # between-cluster predictor mean
  phi_w <- params$phi_w # within-cluster predictor covariance matrix
  phi_b <- params$phi_b # between-cluster predictor covariance matrix

  # get model dimensions
  numXlev1 <- nrow(phi_w)
  numXlev2 <- nrow(phi_b)
  numRaneff <- nrow(cov_b0b1)

  ################## simulate data ##################

  # sample size
  N <- numL2 * numperL2

  # generate level-2 id vector (assumes equal cluster sizes)
  L2id <- seq(1:numL2) %x% matrix(1, nrow = numperL2, ncol = 1)

  # generate within- and between-cluster predictors
  X_w <- rmvnorm(N, mean = rep(0, numXlev1), sigma = phi_w)
  X_b <- rmvnorm(numL2, mean = c(mean_X, mean_Z), sigma = phi_b) %x% matrix(1, nrow = numperL2, ncol = 1)
  X_t <- X_w + X_b

  # generate level-2 random intercept and random slope residuals
  b0b1 <- rmvnorm(numL2, rep(0, numRaneff), cov_b0b1)
  b0b1_stack <- b0b1 %x% matrix(1, nrow = numperL2, ncol = 1)
  if (numRaneff == 1) {
    b0b1_stack <- cbind(b0b1_stack, matrix(0, nrow = nrow(b0b1_stack), ncol = numXlev1))
  }

  # generate level-1 residuals
  e <- rnorm(N, 0, sqrt(var_e))

  # generate y
  Y <- cbind(1, X_w, X_b) %*% betas + rowSums(cbind(1, X_w) * b0b1_stack) + e

  # arrange dataN
  mlmdat <- data.frame(cbind(L2id, Y, X_w, X_b, X_t))
  colnames(mlmdat) <- c("id", "Y", "X1w", "X2w", "X1b", "X2b", "X1t", "X2t")

  return(mlmdat)
}

#################################################################
# perform simulation
#################################################################

# conditions
set.seed(90291)
numL2 <- c(25, 50, 100)
numperL2 <- c(5, 10, 25, 50)
reps <- 1

# get focal model's true values
true_null <- c(params_null$gammas[1:2], params_null$tau, 0, 0, params_null$var_e)
true_alt <- c(params_alt$gammas[1:2], params_alt$tau[1, 1], params_alt$tau[2, 1], params_alt$tau[2, 2], params_alt$var_e)

# iterate
for (l2 in numL2) {
  for (l1 in numperL2) {
    for (r in 1:reps) {
      # simulate one data set
      mlmdata <- gen1data(params = params_ranslope, numL2 = l2, numperL2 = l1)
    }
  }
}
