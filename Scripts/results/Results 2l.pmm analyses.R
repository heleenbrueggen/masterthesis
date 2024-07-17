###########################
# Results 2l.pmm analyses # 
###########################
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
################
# Setting path # 
################
path <- "/Volumes/Heleen 480GB/MBART-MICE files/"
#######################
# Defining parameters #
#######################
combinations <- read_rds(paste(path, "data/combinations.rds", sep = ""))
#############################
# Storing names of datasets #
#############################
names <- read_rds(paste(path, "data/names.rds", sep = ""))
###################
# Loading results #
###################
results_2l.pmm <- list()
for (i in seq_len(nrow(combinations))) {
  results_2l.pmm[[i]] <- read_rds(paste(path, "results/imputed/2l.pmm/analyses_2l.pmm_", names[i], ".rds", sep = ""))
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
  # Obtaining mean of estimates
  mean_estimates <- estimates %>%
    do.call(rbind, .) %>%
    colMeans() %>%
    t() %>%
    as_tibble()
  # Calculating MCSE of bias
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
  # Calculating coverage for all datasets
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
   # Calculating MCSE of coverage
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
  # Calculating confidence interval width for all datasets
  ciw.datasets <- map(estimates, ~ .x %>%
    mutate(ciw = `97.5 %` - `2.5 %`) %>%
    select(ciw) %>%
    t() %>%
    as_tibble() %>%
    rename(beta0j = V1, beta1j = V2, beta2j = V3, beta3j = V4, beta4j = V5, beta5j = V6, beta6j = V7, beta7j = V8, z1 = V9, z2 = V10, `x1:z1` = V11, `x2:z1` = V12, `x3:z2` = V13)) %>%
    list_rbind()
  # Average confidence interval width
  ciw <- colMeans(ciw.datasets) %>%
    t() %>%
    as_tibble()

  return(list(ciw.datasets = ciw.datasets, ciw = ciw))
}
##############
# Evaluation # 
##############
# Bias
bias.datasets_2l.pmm <- list()
for (i in seq_len(nrow(combinations))) {
    # Logging iteration
    cat("Processing iteration:", i, "\n")
    # Bias
    bias.datasets_2l.pmm[[i]] <- bias(results_2l.pmm[[i]])
}
# Extracting relevant information
bias_2l.pmm <- map(bias.datasets_2l.pmm, ~ .x$bias) %>% list_rbind()
mcse_bias_2l.pmm <- map(bias.datasets_2l.pmm, ~ .x$bias.mcse) %>% list_rbind()
# Saving
write_rds(bias.datasets_2l.pmm, file = paste(path, "results/evaluations/bias_datasets_2l.pmm.rds", sep = ""))
write_rds(bias_2l.pmm, file = paste(path, "results/evaluations/bias_2l.pmm.rds", sep = ""))
write_rds(mcse_bias_2l.pmm, file = paste(path, "results/evaluations/mcse_bias_2l.pmm.rds", sep = ""))
# MSE 
mse.datasets_2l.pmm <- list()
for (i in seq_len(nrow(combinations))) {
    # Logging iteration
    cat("Processing iteration:", i, "\n")
    # MSE
    mse.datasets_2l.pmm[[i]] <- mse(results_2l.pmm[[i]])
}
# Extracting relevant information
mse_2l.pmm <- map(mse.datasets_2l.pmm, ~ .x$mse) %>% list_rbind()
mcse_mse_2l.pmm <- map(mse.datasets_2l.pmm, ~ .x$mse.mcse) %>% list_rbind()
# Saving
write_rds(mse.datasets_2l.pmm, file = paste(path, "results/evaluations/mse_datasets_2l.pmm.rds", sep = ""))
write_rds(mse_2l.pmm, file = paste(path, "results/evaluations/mse_2l.pmm.rds", sep = ""))
write_rds(mcse_mse_2l.pmm, file = paste(path, "results/evaluations/mcse_mse_2l.pmm.rds", sep = ""))
# Coverage
coverage.datasets_2l.pmm <- list()
for (i in seq_len(nrow(combinations))) {
    # Logging iteration
    cat("Processing iteration:", i, "\n")
    # Coverage
    coverage.datasets_2l.pmm[[i]] <- coverage(results_2l.pmm[[i]])
}
# Extracting relevant information
coverage_2l.pmm <- map(coverage.datasets_2l.pmm, ~ .x$coverage) %>% list_rbind()
mcse_coverage_2l.pmm <- map(coverage.datasets_2l.pmm, ~ .x$coverage.mcse) %>% list_rbind()
# Saving
write_rds(coverage.datasets_2l.pmm, file = paste(path, "results/evaluations/coverage_datasets_2l.pmm.rds", sep = ""))
write_rds(coverage_2l.pmm, file = paste(path, "results/evaluations/coverage_2l.pmm.rds", sep = ""))
write_rds(mcse_coverage_2l.pmm, file = paste(path, "results/evaluations/mcse_coverage_2l.pmm.rds", sep = ""))
# Confidence interval width
ciw.datasets_2l.pmm <- list()
for (i in seq_len(nrow(combinations))) {
    # Logging iteration
    cat("Processing iteration:", i, "\n")
    # CIW
    ciw.datasets_2l.pmm[[i]] <- ciw(results_2l.pmm[[i]])
}
# Extracting relevant information
ciw_2l.pmm <- map(ciw.datasets_2l.pmm, ~ .x$ciw) %>% list_rbind()
# Saving
write_rds(ciw.datasets_2l.pmm, file = paste(path, "results/evaluations/ciw_datasets_2l.pmm.rds", sep = ""))
write_rds(ciw_2l.pmm, file = paste(path, "results/evaluations/ciw_2l.pmm.rds", sep = ""))