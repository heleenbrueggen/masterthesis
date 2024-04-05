######################################
# Results listwise deletion analyses # 
######################################
#############
# Libraries # 
#############
library(tidyverse)
library(broom.mixed)
################
# Setting seed #
################ 
set.seed(123)
#######################
# Defining parameters #
#######################
# ngroups <- c(30, 50)
# groupsizes <- c(15, 35, 50)
# iccs <- c(.2, .5)
# mar_mcar <- c("mcar")
# miss <- c(0)
# g <- c(.2, .5)
# combinations <- expand.grid(
#   ngroup = ngroups,
#   groupsize = groupsizes,
#   icc = iccs,
#   mar_mcar = mar_mcar,
#   miss = miss,
#   g = g
# )
ngroups <- c(30, 50)
groupsizes <- c(15, 50)
iccs <- c(.5)
mar_mcar <- c("mar", "mcar")
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
#############################
# Storing names of datasets #
#############################
names <- rep(NA, nrow(combinations))
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
###################
# Loading results #
###################
results_ld <- list()
for (i in seq_len(nrow(combinations))) {
  results_ld[[i]] <- read_rds(paste("/Volumes/Heleen\ 480GB/Master\ thesis/results/listwise/analyses_ld_", names[i], ".rds", sep = ""))[1:100]
}
#####################################
# Creating functions for evaluation # 
#####################################
# Bias function
bias <- function(estimated) {
  estimates <- estimated %>%
    map(function(x) {
      fixed <- x %>%
        filter(term == "(Intercept)" | term == "x1" | term == "x2" | term == "x3" | term == "x4" | term == "x5" | term == "x6" | term == "x7" | term == "z1" | term == "z2" | term == "x1:z1" | term == "x2:z1" | term == "x3:z2") %>%
        pull(estimate, term)
      random <- x %>%
        filter(term == "sd__(Intercept)" | term == "sd__x1" | term == "sd__x2" | term == "sd__x3" | term == "sd__Observation") %>%
        pull(estimate, term) %>%
        sapply(function(x) x^2)
      c(fixed, random)
    })

  truth <- tibble(
    beta0j = 10,
    beta1j = .5,
    beta2j = .5,
    beta3j = .5,
    beta4j = .5,
    beta5j = .5,
    beta6j = .5,
    beta7j = .5,
    z1 = .5,
    z2 = .5,
    `x1:z1` = .35,
    `x2:z1` = .35,
    `x3:z2` = .35,
    u0 = 84.74903,
    u1 = 1,
    u2 = 1,
    u3 = 1,
    eij = 25
  )

  bias.datasets <- map(estimates, \(x) (x - truth)) %>%
    list_rbind() %>%
    as_tibble()
  bias <- colMeans(bias.datasets) %>%
    t() %>%
    as_tibble()

  mean_estimates <- estimates %>%
    do.call(rbind, .) %>%
    colMeans() %>%
    t() %>%
    as_tibble()
  mcse <- map(estimates, \(x) (x - mean_estimates)^2) %>%
    list_rbind() %>%
    colSums() %>%
    map_vec(\(x) sqrt(x / (length(estimates) * (length(estimates) - 1)))) %>% 
    t() %>% 
    as_tibble()
  colnames(mcse) <- colnames(bias)

  return(list(bias.datasets = bias.datasets, bias = bias, bias.mcse = mcse))
}
# MSE function
mse <- function(estimated) {
  estimates <- estimated %>%
    map(function(x) {
      fixed <- x %>%
        filter(term == "(Intercept)" | term == "x1" | term == "x2" | term == "x3" | term == "x4" | term == "x5" | term == "x6" | term == "x7" | term == "z1" | term == "z2" | term == "x1:z1" | term == "x2:z1" | term == "x3:z2") %>%
        pull(estimate, term)
      random <- x %>%
        filter(term == "sd__(Intercept)" | term == "sd__x1" | term == "sd__x2" | term == "sd__x3" | term == "sd__Observation") %>%
        pull(estimate, term) %>%
        sapply(function(x) x^2)
      c(fixed, random)
    })

  truth <- tibble(
    beta0j = 10,
    beta1j = .5,
    beta2j = .5,
    beta3j = .5,
    beta4j = .5,
    beta5j = .5,
    beta6j = .5,
    beta7j = .5,
    z1 = .5,
    z2 = .5,
    `x1:z1` = .35,
    `x2:z1` = .35,
    `x3:z2` = .35,
    u0 = 84.74903,
    u1 = 1,
    u2 = 1,
    u3 = 1,
    eij = 25
  )

  mse.datasets <- map(estimates, \(x) (x - truth)^2) %>%
    list_rbind() %>%
    as_tibble()
  mse <- colMeans(mse.datasets) %>%
    t() %>%
    as_tibble()

  mcse <- map(estimates, \(x) (x - truth)^2) %>%
    map(., \(x) (x - mse)^2) %>%
    list_rbind() %>%
    colSums() %>%
    map_vec(\(x) sqrt(x / (length(estimates) * (length(estimates) - 1)))) %>% 
    t() %>%
    as_tibble()

  return(list(mse.datasets = mse.datasets, mse = mse, mse.mcse = mcse))
}
# Coverage function
coverage <- function(estimated) {
  estimates <- estimated %>%
    map(~ .x %>%
      filter(term == "(Intercept)" | term == "x1" | term == "x2" | term == "x3" | term == "x4" | term == "x5" | term == "x6" | term == "x7" | term == "z1" | term == "z2" | term == "x1:z1" | term == "x2:z1" | term == "x3:z2") %>%
      select(conf.low, conf.high))

  truth <- tibble(
    beta0j = 10,
    beta1j = .5,
    beta2j = .5,
    beta3j = .5,
    beta4j = .5,
    beta5j = .5,
    beta6j = .5,
    beta7j = .5,
    z1 = .5,
    z2 = .5,
    `x1:z1` = .35,
    `x2:z1` = .35,
    `x3:z2` = .35
  ) %>% 
    t() %>% 
    as_tibble() %>% 
    rename(value = V1)

  combined <- map(estimates, \(x) cbind(x, truth))

  coverage.datasets <- map(combined, \(x) x %>%
    mutate(coverage = value > conf.low & value < conf.high) %>%
    select(coverage) %>%
    mutate(coverage = as.numeric(coverage)) %>%
    t() %>%
    as_tibble() %>%
    rename(beta0j = V1, beta1j = V2, beta2j = V3, beta3j = V4, beta4j = V5, beta5j = V6, beta6j = V7, beta7j = V8, z1 = V9, z2 = V10, `x1:z1` = V11, `x2:z1` = V12, `x3:z2` = V13)) %>%
    list_rbind()

  coverage <- colMeans(coverage.datasets) %>%
    t() %>%
    as_tibble()

  mcse <- sqrt(coverage * (1 - coverage) / length(estimates)) %>%
    as_tibble()

  return(list(coverage.datasets = coverage.datasets, coverage = coverage, coverage.mcse = mcse))
}
# Confidence interval width function
ciw <- function(estimated) { 
  estimates <- estimated %>%
    map(~ .x %>%
      filter(term == "(Intercept)" | term == "x1" | term == "x2" | term == "x3" | term == "x4" | term == "x5" | term == "x6" | term == "x7" | term == "z1" | term == "z2" | term == "x1:z1" | term == "x2:z1" | term == "x3:z2") %>%
      select(conf.low, conf.high))

    ciw.datasets <- map(estimates, ~ .x %>%
      mutate(ciw = conf.high - conf.low) %>%
      select(ciw) %>%
      t() %>%
      as_tibble() %>%
      rename(beta0j = V1, beta1j = V2, beta2j = V3, beta3j = V4, beta4j = V5, beta5j = V6, beta6j = V7, beta7j = V8, z1 = V9, z2 = V10, `x1:z1` = V11, `x2:z1` = V12, `x3:z2` = V13)) %>% 
      list_rbind()

    ciw <- colMeans(ciw.datasets) %>%
      t() %>%
      as_tibble()
    
    return(list(ciw.datasets = ciw.datasets, ciw = ciw))
}
##############
# Evaluation # 
##############
# Bias
bias.datasets_ld <- list()
for (i in seq_len(nrow(combinations))) {
    # Logging iteration
    cat("Processing iteration:", i, "\n")
    # Bias
    bias.datasets_ld[[i]] <- bias(results_ld[[i]])
}
bias_ld <- map(bias.datasets_ld, ~ .x$bias) %>% list_rbind()
mcse_bias_ld <- map(bias.datasets_ld, ~ .x$bias.mcse) %>% list_rbind()
# Saving
write_rds(bias.datasets_ld, file = "/Volumes/Heleen\ 480GB/Master\ thesis/results/evaluations/bias_datasets_ld.rds")
write_rds(bias_ld, file = "/Volumes/Heleen\ 480GB/Master\ thesis/results/evaluations/bias_ld.rds")
write_rds(mcse_bias_ld, file = "/Volumes/Heleen\ 480GB/Master\ thesis/results/evaluations/mcse_bias_ld.rds")
# MSE 
mse.datasets_ld <- list()
for (i in seq_len(nrow(combinations))) {
    # Logging iteration
    cat("Processing iteration:", i, "\n")
    # MSE
    mse.datasets_ld[[i]] <- mse(results_ld[[i]])
}
mse_ld <- map(mse.datasets_ld, ~ .x$mse) %>% list_rbind()
mcse_mse_ld <- map(mse.datasets_ld, ~ .x$mse.mcse) %>% list_rbind()
# Saving
write_rds(mse.datasets_ld, file = "/Volumes/Heleen\ 480GB/Master\ thesis/results/evaluations/mse_datasets_ld.rds")
write_rds(mse_ld, file = "/Volumes/Heleen\ 480GB/Master\ thesis/results/evaluations/mse_ld.rds")
write_rds(mcse_mse_ld, file = "/Volumes/Heleen\ 480GB/Master\ thesis/results/evaluations/mcse_mse_ld.rds")
# Coverage
coverage.datasets_ld <- list()
for (i in seq_len(nrow(combinations))) {
    # Logging iteration
    cat("Processing iteration:", i, "\n")
    # Coverage
    coverage.datasets_ld[[i]] <- coverage(results_ld[[i]])
}
coverage_ld <- map(coverage.datasets_ld, ~ .x$coverage) %>% list_rbind()
mcse_coverage_ld <- map(coverage.datasets_ld, ~ .x$coverage.mcse) %>% list_rbind()
# Saving
write_rds(coverage.datasets_ld, file = "/Volumes/Heleen\ 480GB/Master\ thesis/results/evaluations/coverage_datasets_ld.rds")
write_rds(coverage_ld, file = "/Volumes/Heleen\ 480GB/Master\ thesis/results/evaluations/coverage_ld.rds")
write_rds(mcse_coverage_ld, file = "/Volumes/Heleen\ 480GB/Master\ thesis/results/evaluations/mcse_coverage_ld.rds")
# CIW 
ciw.datasets_ld <- list()
for (i in seq_len(nrow(combinations))) {
    # Logging iteration
    cat("Processing iteration:", i, "\n")
    # CIW
    ciw.datasets_ld[[i]] <- ciw(results_ld[[i]])
}
ciw_ld <- map(ciw.datasets_ld, ~ .x$ciw) %>% list_rbind()
# Saving
write_rds(ciw.datasets_ld, file = "/Volumes/Heleen\ 480GB/Master\ thesis/results/evaluations/ciw_datasets_ld.rds")
write_rds(ciw_ld, file = "/Volumes/Heleen\ 480GB/Master\ thesis/results/evaluations/ciw_ld.rds")