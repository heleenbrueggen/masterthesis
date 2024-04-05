########################
# Results pmm analyses # 
########################
#############
# Libraries # 
#############
library(tidyverse)
library(broom.mixed)
library(mitml)
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
# mar_mcar <- c("mar", "mcar")
# miss <- c(25, 50)
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
results_pmm <- list()
for (i in seq_len(nrow(combinations))) {
  results_pmm[[i]] <- read_rds(paste("/Volumes/Heleen\ 480GB/Master\ thesis/results/imputed/pmm/analyses_pmm_", names[i], ".rds", sep = ""))[1:100]
}
#####################################
# Creating functions for evaluation # 
#####################################
# Bias function
bias <- function(estimated) {
  estimates <- estimated %>%
    map(., \(x) c(
      x[2]$estimates[1:13, 1],
      x[3]$extra.pars[c(1, 2, 3, 4, 11), ]
    ))

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
    map(., \(x) c(
      x[2]$estimates[1:13, 1],
      x[3]$extra.pars[c(1, 2, 3, 4, 11), ]
    ))

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
    map(~ .x %>% mitml::confint.mitml.testEstimates())

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

  coverage.datasets <- map(combined, ~ .x %>%
    mutate(coverage = value > `2.5 %` & value < `97.5 %`) %>%
    select(coverage) %>%
    mutate(coverage = as.numeric(coverage)) %>%
    t() %>%
    as_tibble() %>%
    rename(beta0j = `(Intercept)`, beta1j = x1, beta2j = x2, beta3j = x3, beta4j = x4, beta5j = x5, beta6j = x6, beta7j = x7, z1 = z1, z2 = z2, `x1:z1` = `x1:z1`, `x2:z1` = `x2:z1`, `x3:z2` = `x3:z2`)) %>%
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
    map(~ .x %>% mitml::confint.mitml.testEstimates() %>% 
    as_tibble())

    ciw.datasets <- map(estimates, ~ .x %>%
      mutate(ciw = `97.5 %` - `2.5 %`) %>%
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
bias.datasets_pmm <- list()
for (i in seq_len(nrow(combinations))) {
    # Logging iteration
    cat("Processing iteration:", i, "\n")
    # Bias
    bias.datasets_pmm[[i]] <- bias(results_pmm[[i]])
}
bias_pmm <- map(bias.datasets_pmm, ~ .x$bias) %>% list_rbind()
mcse_bias_pmm <- map(bias.datasets_pmm, ~ .x$bias.mcse) %>% list_rbind()
# Saving
write_rds(bias.datasets_pmm, file = "/Volumes/Heleen\ 480GB/Master\ thesis/results/evaluations/bias_datasets_pmm.rds")
write_rds(bias_pmm, file = "/Volumes/Heleen\ 480GB/Master\ thesis/results/evaluations/bias_pmm.rds")
write_rds(mcse_bias_pmm, file = "/Volumes/Heleen\ 480GB/Master\ thesis/results/evaluations/mcse_bias_pmm.rds")
# MSE 
mse.datasets_pmm <- list()
for (i in seq_len(nrow(combinations))) {
    # Logging iteration
    cat("Processing iteration:", i, "\n")
    # MSE
    mse.datasets_pmm[[i]] <- mse(results_pmm[[i]])
}
mse_pmm <- map(mse.datasets_pmm, ~ .x$mse) %>% list_rbind()
mcse_mse_pmm <- map(mse.datasets_pmm, ~ .x$mse.mcse) %>% list_rbind()
# Saving
write_rds(mse.datasets_pmm, file = "/Volumes/Heleen\ 480GB/Master\ thesis/results/evaluations/mse_datasets_pmm.rds")
write_rds(mse_pmm, file = "/Volumes/Heleen\ 480GB/Master\ thesis/results/evaluations/mse_pmm.rds")
write_rds(mcse_mse_pmm, file = "/Volumes/Heleen\ 480GB/Master\ thesis/results/evaluations/mcse_mse_pmm.rds")
# Coverage
coverage.datasets_pmm <- list()
for (i in seq_len(nrow(combinations))) {
    # Logging iteration
    cat("Processing iteration:", i, "\n")
    # Coverage
    coverage.datasets_pmm[[i]] <- coverage(results_pmm[[i]])
}
coverage_pmm <- map(coverage.datasets_pmm, ~ .x$coverage) %>% list_rbind()
mcse_coverage_pmm <- map(coverage.datasets_pmm, ~ .x$coverage.mcse) %>% list_rbind()
# Saving
write_rds(coverage.datasets_pmm, file = "/Volumes/Heleen\ 480GB/Master\ thesis/results/evaluations/coverage_datasets_pmm.rds")
write_rds(coverage_pmm, file = "/Volumes/Heleen\ 480GB/Master\ thesis/results/evaluations/coverage_pmm.rds")
write_rds(mcse_coverage_pmm, file = "/Volumes/Heleen\ 480GB/Master\ thesis/results/evaluations/mcse_coverage_pmm.rds")
# CIW
ciw.datasets_pmm <- list()
for (i in seq_len(nrow(combinations))) {
    # Logging iteration
    cat("Processing iteration:", i, "\n")
    # CIW
    ciw.datasets_pmm[[i]] <- ciw(results_pmm[[i]])
}
ciw_pmm <- map(ciw.datasets_pmm, ~ .x$ciw) %>% list_rbind()
# Saving
write_rds(ciw.datasets_pmm, file = "/Volumes/Heleen\ 480GB/Master\ thesis/results/evaluations/ciw_datasets_pmm.rds")
write_rds(ciw_pmm, file = "/Volumes/Heleen\ 480GB/Master\ thesis/results/evaluations/ciw_pmm.rds")