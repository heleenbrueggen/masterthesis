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
results_2l.pmm <- list()
for (i in seq_len(nrow(combinations))) {
  results_2l.pmm[[i]] <- read_rds(paste("/Volumes/Heleen\ 480GB/Master\ thesis/results/imputed/2l.pmm/analyses_2l.pmm_", names[i], ".rds", sep = ""))[1:100]
}
#####################################
# Creating functions for evaluation # 
#####################################
# Bias function
bias <- function(estimated, true) {
  estimates <- estimated %>%
    map(., \(x) c(
      x[2]$estimates[1:13, 1],
      x[3]$extra.pars[c(1, 2, 3, 4, 11), ]
    ))

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

  bias.datasets <- map2(estimates, truth, ~ .x - .y) %>%
    list_rbind() %>%
    as_tibble()
  bias <- colMeans(bias.datasets) %>%
    t() %>%
    as_tibble()

  return(list(bias.datasets = bias.datasets, bias = bias))
}
# MSE function
mse <- function(estimated, true) {
  estimates <- estimated %>%
    map(., \(x) c(
      x[2]$estimates[1:13, 1],
      x[3]$extra.pars[c(1, 2, 3, 4, 11), ]
    ))

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
        across(c(u0, u1, u2, u3, eij), sd)
      ))

  mse.datasets <- map2(estimates, truth, ~ (.x - .y) ^ 2) %>% list_rbind() %>% as_tibble()
  mse <- colMeans(mse.datasets) %>% t() %>% as_tibble()

  return(list(mse.datasets = mse.datasets, mse = mse))
}
coverage <- function(estimated, true) {
  estimates <- estimated %>%
    map(~ .x %>% mitml::confint.mitml.testEstimates())

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
    mutate(coverage = value > `2.5 %` & value < `97.5 %`) %>% 
    select(coverage) %>% 
    mutate(coverage = as.numeric(coverage)) %>% 
    t() %>% 
    as_tibble() %>% 
    rename(beta0j = `(Intercept)`, beta1j = x1, beta2j = x2, beta3j = x3, beta4j = x4, beta5j = x5, beta6j = x6, beta7j = x7, z1 = z1, z2 = z2, `x1:z1` = `x1:z1`, `x2:z1` = `x2:z1`, `x3:z2` = `x3:z2`)) %>% 
    list_rbind()

  coverage <- colMeans(coverage.datasets) %>% t() %>% as_tibble()

  return(list(coverage.datasets = coverage.datasets, coverage = coverage))
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
    bias.datasets_2l.pmm[[i]] <- bias(results_2l.pmm[[i]], read_rds(paste("/Volumes/Heleen\ 480GB/Master\ thesis/data/complete/simdata_", names[i], ".rds", sep = ""))[1:100])
}
bias_2l.pmm <- map(bias.datasets_2l.pmm, ~ .x$bias) %>% list_rbind()
# Saving
write_rds(bias.datasets_2l.pmm, file = "/Volumes/Heleen\ 480GB/Master\ thesis/results/evaluations/bias_datasets_2l.pmm.rds")
write_rds(bias_2l.pmm, file = "/Volumes/Heleen\ 480GB/Master\ thesis/results/evaluations/bias_2l.pmm.rds")
# MSE 
mse.datasets_2l.pmm <- list()
for (i in seq_len(nrow(combinations))) {
    # Logging iteration
    cat("Processing iteration:", i, "\n")
    # MSE
    mse.datasets_2l.pmm[[i]] <- mse(results_2l.pmm[[i]], read_rds(paste("/Volumes/Heleen\ 480GB/Master\ thesis/data/complete/simdata_", names[i], ".rds", sep = ""))[1:100])
}
mse_2l.pmm <- map(mse.datasets_2l.pmm, ~ .x$mse) %>% list_rbind()
# Saving
write_rds(mse.datasets_2l.pmm, file = "/Volumes/Heleen\ 480GB/Master\ thesis/results/evaluations/mse_datasets_2l.pmm.rds")
write_rds(mse_2l.pmm, file = "/Volumes/Heleen\ 480GB/Master\ thesis/results/evaluations/mse_2l.pmm.rds")
# Coverage
coverage.datasets_2l.pmm <- list()
for (i in seq_len(nrow(combinations))) {
    # Logging iteration
    cat("Processing iteration:", i, "\n")
    # Coverage
    coverage.datasets_2l.pmm[[i]] <- coverage(results_2l.pmm[[i]], read_rds(paste("/Volumes/Heleen\ 480GB/Master\ thesis/data/complete/simdata_", names[i], ".rds", sep = ""))[1:100])
}
coverage_2l.pmm <- map(coverage.datasets_2l.pmm, ~ .x$coverage) %>% list_rbind()
# Saving
write_rds(coverage.datasets_2l.pmm, file = "/Volumes/Heleen\ 480GB/Master\ thesis/results/evaluations/coverage_datasets_2l.pmm.rds")
write_rds(coverage_2l.pmm, file = "/Volumes/Heleen\ 480GB/Master\ thesis/results/evaluations/coverage_2l.pmm.rds")
