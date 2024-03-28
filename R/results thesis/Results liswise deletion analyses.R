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
mar_mcar <- c("mar")
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
bias <- function(estimated, true) {
  estimates <- estimated %>% 
  map(function(x) { 
    fixed <- x %>% filter(term == "(Intercept)" | term == "x1" | term == "x2" | term == "x3" | term == "x4" | term == "x5" | term == "x6" | term == "x7" | term == "z1" | term == "z2" | term == "x1:z1" | term == "x2:z1" | term == "x3:z2") %>% pull(estimate, term)
    random <- x %>% filter(term == "sd__(Intercept)" | term == "sd__x1" | term == "sd__x2" | term == "sd__x3" | term == "sd__Observation") %>% pull(estimate, term) %>% sapply(function(x) x^2)
    c(fixed, random)
  })

  truth <- true %>%
    map(~ .x %>%
      select(beta0j, beta1j, beta2j, beta3j, beta4j, beta5j, beta6j, beta7j, u0, u1, u2, u3, eij) %>%
      mutate(
        z1 = .5,
        z2 = .5,
        `x1:z1` = .35,
        `x2:z1` = .35,
        `x3:z2` = .35
      ) %>%
      relocate(c(z1, z2, `x1:z1`, `x2:z1`, `x3:z2`), .before = "u0") %>%
      summarise(
        across(c(beta0j, beta1j, beta2j, beta3j, beta4j, beta5j, beta6j, beta7j, z1, z2, `x1:z1`, `x2:z1`, `x3:z2`), mean),
        across(c(u0, u1, u2, u3, eij), var)
      ))

  bias.datasets <- map2(estimates, truth, ~ .x - .y) %>% list_rbind() %>% as_tibble()
  bias <- colMeans(bias.datasets) %>% t() %>% as_tibble()

  return(list(bias.datasets = bias.datasets, bias = bias))
}
# MSE function
mse <- function(estimated, true) {
  estimates <- estimated %>% 
  map(function(x) { 
    fixed <- x %>% filter(term == "(Intercept)" | term == "x1" | term == "x2" | term == "x3" | term == "x4" | term == "x5" | term == "x6" | term == "x7" | term == "z1" | term == "z2" | term == "x1:z1" | term == "x2:z1" | term == "x3:z2") %>% pull(estimate, term)
    random <- x %>% filter(term == "sd__(Intercept)" | term == "sd__x1" | term == "sd__x2" | term == "sd__x3" | term == "sd__Observation") %>% pull(estimate, term) %>% sapply(function(x) x^2)
    c(fixed, random)
  })

  truth <- true %>%
    map(~ .x %>%
      select(beta0j, beta1j, beta2j, beta3j, beta4j, beta5j, beta6j, beta7j, u0, u1, u2, u3, eij) %>%
      mutate(
        z1 = .5,
        z2 = .5,
        `x1:z1` = .35,
        `x2:z1` = .35,
        `x3:z2` = .35
      ) %>%
      relocate(c(z1, z2, `x1:z1`, `x2:z1`, `x3:z2`), .before = "u0") %>%
      summarise(
        across(c(beta0j, beta1j, beta2j, beta3j, beta4j, beta5j, beta6j, beta7j, z1, z2, `x1:z1`, `x2:z1`, `x3:z2`), mean),
        across(c(u0, u1, u2, u3, eij), var)
      ))

  mse.datasets <- map2(estimates, truth, ~ (.x - .y) ^ 2) %>% list_rbind() %>% as_tibble()
  mse <- colMeans(mse.datasets) %>% t() %>% as_tibble()

  return(list(mse.datasets = mse.datasets, mse = mse))
}
# Coverage function
coverage <- function(estimated, true) {
  estimates <- estimated %>%
    map(~ .x %>%
      filter(term == "(Intercept)" | term == "x1" | term == "x2" | term == "x3" | term == "x4" | term == "x5" | term == "x6" | term == "x7" | term == "z1" | term == "z2" | term == "x1:z1" | term == "x2:z1" | term == "x3:z2") %>%
      select(conf.low, conf.high))

  truth <- true %>%
    map(~ .x %>%
      select(beta0j, beta1j, beta2j, beta3j, beta4j, beta5j, beta6j, beta7j) %>%
      mutate(
        z1 = .5, 
        z2 = .5, 
        `x1:z1` = .35, 
        `x2:z1` = .35, 
        `x3:z2` = .35) %>%
      summarise(
        across(c(beta0j, beta1j, beta2j, beta3j, beta4j, beta5j, beta6j, beta7j, z1, z2, `x1:z1`, `x2:z1`, `x3:z2`), mean)) %>% 
        pivot_longer(everything(), names_to = "term", values_to = "value"))

  combined <- map2(estimates, truth, ~ .y %>% cbind(.x))

  coverage.datasets <- map(combined, ~ .x %>%
    mutate(coverage = value > conf.low & value < conf.high) %>% 
    select(coverage) %>% 
    mutate(coverage = as.numeric(coverage)) %>% 
    t() %>% 
    as_tibble() %>% 
    rename(beta0j = V1, beta1j = V2, beta2j = V3, beta3j = V4, beta4j = V5, beta5j = V6, beta6j = V7, beta7j = V8, z1 = V9, z2 = V10, `x1:z1` = V11, `x2:z1` = V12, `x3:z2` = V13)) %>% 
    list_rbind()

  coverage <- colMeans(coverage.datasets) %>% t() %>% as_tibble()

  return(list(coverage.datasets = coverage.datasets, coverage = coverage))
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
    bias.datasets_ld[[i]] <- bias(results_ld[[i]], read_rds(paste("/Volumes/Heleen\ 480GB/Master\ thesis/data/complete/simdata_", names[i], ".rds", sep = ""))[1:100])
}
bias_ld <- map(bias.datasets_ld, ~ .x$bias) %>% list_rbind()
# Saving
write_rds(bias.datasets_ld, file = "/Volumes/Heleen\ 480GB/Master\ thesis/results/evaluations/bias_datasets_ld.rds")
write_rds(bias_ld, file = "/Volumes/Heleen\ 480GB/Master\ thesis/results/evaluations/bias_ld.rds")
# MSE 
mse.datasets_ld <- list()
for (i in seq_len(nrow(combinations))) {
    # Logging iteration
    cat("Processing iteration:", i, "\n")
    # MSE
    mse.datasets_ld[[i]] <- mse(results_ld[[i]], read_rds(paste("/Volumes/Heleen\ 480GB/Master\ thesis/data/complete/simdata_", names[i], ".rds", sep = ""))[1:100])
}
mse_ld <- map(mse.datasets_ld, ~ .x$mse) %>% list_rbind()
# Saving
write_rds(mse.datasets_ld, file = "/Volumes/Heleen\ 480GB/Master\ thesis/results/evaluations/mse_datasets_ld.rds")
write_rds(mse_ld, file = "/Volumes/Heleen\ 480GB/Master\ thesis/results/evaluations/mse_ld.rds")
# Coverage
coverage.datasets_ld <- list()
for (i in seq_len(nrow(combinations))) {
    # Logging iteration
    cat("Processing iteration:", i, "\n")
    # Coverage
    coverage.datasets_ld[[i]] <- coverage(results_ld[[i]], read_rds(paste("/Volumes/Heleen\ 480GB/Master\ thesis/data/complete/simdata_", names[i], ".rds", sep = ""))[1:100])
}
coverage_ld <- map(coverage.datasets_ld, ~ .x$coverage) %>% list_rbind()
# Saving
write_rds(coverage.datasets_ld, file = "/Volumes/Heleen\ 480GB/Master\ thesis/results/evaluations/coverage_datasets_ld.rds")
write_rds(coverage_ld, file = "/Volumes/Heleen\ 480GB/Master\ thesis/results/evaluations/coverage_ld.rds")
