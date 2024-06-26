########################
# Results stan4bart analyses # 
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
ngroups <- c(30, 50) # Number of groups 
groupsizes <- c(15, 50) # Group sizes 
iccs <- c(.5) # Intraclass correlation coefficient
mar_mcar <- c("mar", "mcar") # Missing data mechanism
miss <- c(50) # Percentage of missing data
g <- c(.5) # Within-group effect size
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
results_stan4bart <- list()
for (i in seq_len(nrow(combinations))) {
  results_stan4bart[[i]] <- read_rds(paste("results/imputed/stan4bart/analyses_stan4bart_", names[i], ".rds", sep = ""))
}
#####################################
# Creating functions for evaluation # 
#####################################
# Bias function
bias <- function(estimated) {
  # Extracting estimates
  estimates <- estimated %>%
    map(., \(x) c(
      x[2]$estimates[1:13, 1],
      x[3]$extra.pars[c(1, 2, 3, 4, 11), ]
    ))
  # Defininf truth
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
  # Calculating bias for all datasets
  bias.datasets <- map(estimates, \(x) (x - truth)) %>%
    list_rbind() %>%
    as_tibble()
  # Average bias
  bias <- colMeans(bias.datasets) %>%
    t() %>%
    as_tibble()
  # Mean of estimates
  mean_estimates <- estimates %>%
    do.call(rbind, .) %>%
    colMeans() %>%
    t() %>%
    as_tibble()
  # MCSE of bias
  mcse <- map(estimates, \(x) (x - mean_estimates)^2) %>%
    list_rbind() %>%
    colSums() %>%
    map_vec(\(x) sqrt(x / (length(estimates) * (length(estimates) - 1)))) %>%
    t() %>% 
    as_tibble()
  colnames(mcse) <- colnames(bias)

  return(list(bias.datasets = bias.datasets, bias = bias, bias.mcse = mcse))
}
# Coverage function
coverage <- function(estimated) {
  # Extracting estimates
  estimates <- estimated %>%
    map(~ .x %>% mitml::confint.mitml.testEstimates())
  # Defining truth
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
  # Combining estimates and truth
  combined <- map(estimates, \(x) cbind(x, truth))
  # Coverage of all data sets
  coverage.datasets <- map(combined, ~ .x %>%
    mutate(coverage = value > `2.5 %` & value < `97.5 %`) %>%
    select(coverage) %>%
    mutate(coverage = as.numeric(coverage)) %>%
    t() %>%
    as_tibble() %>%
    rename(beta0j = `(Intercept)`, beta1j = x1, beta2j = x2, beta3j = x3, beta4j = x4, beta5j = x5, beta6j = x6, beta7j = x7, z1 = z1, z2 = z2, `x1:z1` = `x1:z1`, `x2:z1` = `x2:z1`, `x3:z2` = `x3:z2`)) %>%
    list_rbind()
  # Average coverage
  coverage <- colMeans(coverage.datasets) %>%
    t() %>%
    as_tibble()
  # MCSE of coverage
  mcse <- sqrt(coverage * (1 - coverage) / length(estimates)) %>%
    as_tibble()

  return(list(coverage.datasets = coverage.datasets, coverage = coverage, coverage.mcse = mcse))
}
# Confidence interval width function
ciw <- function(estimated) {
  # Extracting estimates
  estimates <- estimated %>%
    map(~ .x %>%
      mitml::confint.mitml.testEstimates() %>%
      as_tibble())
  # Calculating CIW for all datasets
  ciw.datasets <- map(estimates, ~ .x %>%
    mutate(ciw = `97.5 %` - `2.5 %`) %>%
    select(ciw) %>%
    t() %>%
    as_tibble() %>%
    rename(beta0j = V1, beta1j = V2, beta2j = V3, beta3j = V4, beta4j = V5, beta5j = V6, beta6j = V7, beta7j = V8, z1 = V9, z2 = V10, `x1:z1` = V11, `x2:z1` = V12, `x3:z2` = V13)) %>%
    list_rbind()
  # Average CIW
  ciw <- colMeans(ciw.datasets) %>%
    t() %>%
    as_tibble()

  return(list(ciw.datasets = ciw.datasets, ciw = ciw))
}
##############
# Evaluation # 
##############
# Bias
bias.datasets_stan4bart <- list()
for (i in seq_len(nrow(combinations))) {
    # Logging iteration
    cat("Processing iteration:", i, "\n")
    # Bias
    bias.datasets_stan4bart[[i]] <- bias(results_stan4bart[[i]])
}
# Extracting relevant information
bias_stan4bart <- map(bias.datasets_stan4bart, ~ .x$bias) %>% list_rbind()
mcse_bias_stan4bart <- map(bias.datasets_stan4bart, ~ .x$bias.mcse) %>% list_rbind()
# Saving
write_rds(bias.datasets_stan4bart, file = "results/evaluations/bias_datasets_stan4bart.rds")
write_rds(bias_stan4bart, file = "results/evaluations/bias_stan4bart.rds")
write_rds(mcse_bias_stan4bart, file = "results/evaluations/mcse_bias_stan4bart.rds")
# MSE 
mse.datasets_stan4bart <- list()
for (i in seq_len(nrow(combinations))) {
    # Logging iteration
    cat("Processing iteration:", i, "\n")
    # MSE
    mse.datasets_stan4bart[[i]] <- mse(results_stan4bart[[i]])
}
# Extracting relevant information
mse_stan4bart <- map(mse.datasets_stan4bart, ~ .x$mse) %>% list_rbind()
mcse_mse_stan4bart <- map(mse.datasets_stan4bart, ~ .x$mse.mcse) %>% list_rbind()
# Saving
write_rds(mse.datasets_stan4bart, file = "results/evaluations/mse_datasets_stan4bart.rds")
write_rds(mse_stan4bart, file = "results/evaluations/mse_stan4bart.rds") 
write_rds(mcse_mse_stan4bart, file = "results/evaluations/mcse_mse_stan4bart.rds")
# Coverage
coverage.datasets_stan4bart <- list()
for (i in seq_len(nrow(combinations))) {
    # Logging iteration
    cat("Processing iteration:", i, "\n")
    # Coverage
    coverage.datasets_stan4bart[[i]] <- coverage(results_stan4bart[[i]])
}
# Extracting relevant information
coverage_stan4bart <- map(coverage.datasets_stan4bart, ~ .x$coverage) %>% list_rbind()
mcse_coverage_stan4bart <- map(coverage.datasets_stan4bart, ~ .x$coverage.mcse) %>% list_rbind()
# Saving
write_rds(coverage.datasets_stan4bart, file = "results/evaluations/coverage_datasets_stan4bart.rds")
write_rds(coverage_stan4bart, file = "results/evaluations/coverage_stan4bart.rds")
write_rds(mcse_coverage_stan4bart, file = "results/evaluations/mcse_coverage_stan4bart.rds")
# CIW
ciw.datasets_stan4bart <- list()
for (i in seq_len(nrow(combinations))) {
    # Logging iteration
    cat("Processing iteration:", i, "\n")
    # CIW
    ciw.datasets_stan4bart[[i]] <- ciw(results_stan4bart[[i]])
}
# Extracting relevant information
ciw_stan4bart <- map(ciw.datasets_stan4bart, ~ .x$ciw) %>% list_rbind()
# Saving
write_rds(ciw.datasets_stan4bart, file = "results/evaluations/ciw_datasets_stan4bart.rds")
write_rds(ciw_stan4bart, file = "results/evaluations/ciw_stan4bart.rds")