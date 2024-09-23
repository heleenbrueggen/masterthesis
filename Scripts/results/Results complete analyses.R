#############################
# Results complete analyses # 
#############################
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
combinations <- read_rds(paste(path, "data/nomissing/combinations.rds", sep = ""))
#############################
# Storing names of datasets #
#############################
names <- read_rds(paste(path, "data/nomissing/names.rds", sep = ""))
###################
# Loading results #
###################
results_nomiss <- list()
for (i in seq_len(nrow(combinations))) {
  results_nomiss[[i]] <- read_rds(paste(path, "results/nomissing/analyses_nomiss_", names[i], ".rds", sep = ""))
}
#####################################
# Creating functions for evaluation # 
#####################################
# Bias function
bias <- function(estimated) {
  # Extracting estimates
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
  # Bias for all data sets
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
  # Calculating MCSE
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
    map(~ .x %>%
      filter(term == "(Intercept)" | term == "x1" | term == "x2" | term == "x3" | term == "x4" | term == "x5" | term == "x6" | term == "x7" | term == "z1" | term == "z2" | term == "x1:z1" | term == "x2:z1" | term == "x3:z2") %>%
      select(conf.low, conf.high))
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
  coverage.datasets <- map(combined, \(x) x %>%
    mutate(coverage = value > conf.low & value < conf.high) %>%
    select(coverage) %>%
    mutate(coverage = as.numeric(coverage)) %>%
    t() %>%
    as_tibble() %>%
    rename(beta0j = V1, beta1j = V2, beta2j = V3, beta3j = V4, beta4j = V5, beta5j = V6, beta6j = V7, beta7j = V8, z1 = V9, z2 = V10, `x1:z1` = V11, `x2:z1` = V12, `x3:z2` = V13)) %>%
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
      filter(term == "(Intercept)" | term == "x1" | term == "x2" | term == "x3" | term == "x4" | term == "x5" | term == "x6" | term == "x7" | term == "z1" | term == "z2" | term == "x1:z1" | term == "x2:z1" | term == "x3:z2") %>%
      select(conf.low, conf.high))
  # Calculating CIW for all datasets
  ciw.datasets <- map(estimates, ~ .x %>%
    mutate(ciw = conf.high - conf.low) %>%
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
bias.datasets_nomiss <- list()
for (i in seq_len(nrow(combinations))) {
    # Logging iteration
    cat("Processing iteration:", i, "\n")
    # Bias
    bias.datasets_nomiss[[i]] <- bias(results_nomiss[[i]])
}
# Extracting relevant information
bias_nomiss <- map(bias.datasets_nomiss, ~ .x$bias) %>% list_rbind()
mcse_bias_nomiss <- map(bias.datasets_nomiss, ~ .x$bias.mcse) %>% list_rbind()
# Saving
write_rds(bias.datasets_nomiss, file = paste(path, "results/evaluations/bias_datasets_nomiss.rds", sep = ""))
write_rds(bias_nomiss, file = paste(path, "results/evaluations/bias_nomiss.rds", sep = ""))
write_rds(mcse_bias_nomiss, file = paste(path, "results/evaluations/mcse_bias_nomiss.rds", sep = ""))
# Coverage
coverage.datasets_nomiss <- list()
for (i in seq_len(nrow(combinations))) {
    # Logging iteration
    cat("Processing iteration:", i, "\n")
    # Coverage
    coverage.datasets_nomiss[[i]] <- coverage(results_nomiss[[i]])
}
# Extracting relevant information
coverage_nomiss <- map(coverage.datasets_nomiss, ~ .x$coverage) %>% list_rbind()
mcse_coverage_nomiss <- map(coverage.datasets_nomiss, ~ .x$coverage.mcse) %>% list_rbind()
# Saving
write_rds(coverage.datasets_nomiss, file = paste(path, "results/evaluations/coverage_datasets_nomiss.rds", sep = ""))
write_rds(coverage_nomiss, file = paste(path, "results/evaluations/coverage_nomiss.rds", sep = ""))
write_rds(mcse_coverage_nomiss, file = paste(path, "results/evaluations/mcse_coverage_nomiss.rds", sep = ""))
# CIW
ciw.datasets_nomiss <- list()
for (i in seq_len(nrow(combinations))) {
    # Logging iteration
    cat("Processing iteration:", i, "\n")
    # CIW
    ciw.datasets_nomiss[[i]] <- ciw(results_nomiss[[i]])
}
# Extracting relevant information
ciw_nomiss <- map(ciw.datasets_nomiss, ~ .x$ciw) %>% list_rbind()
# Saving
write_rds(ciw.datasets_nomiss, file = paste(path, "results/evaluations/ciw_datasets_nomiss.rds", sep = ""))
write_rds(ciw_nomiss, file = paste(path, "results/evaluations/ciw_nomiss.rds", sep = ""))