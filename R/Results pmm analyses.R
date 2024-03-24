########################
# Results pmm analyses # 
########################
#############
# Libraries # 
#############
library(ggplot2)
library(tidyverse)
library(broom.mixed)
library(lme4)
################
# Setting seed #
################ 
set.seed(123)
#######################
# Defining parameters #
#######################
ngroups <- c(30, 50)
groupsizes <- c(15, 35, 50)
iccs <- c(.2, .5)
mar_mcar <- c("mar", "mcar")
miss <- c(25, 50)
g <- c(.2, .5)
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
  results_pmm[[i]] <- map(read_rds(paste("/Volumes/Heleen\ 480GB/Master\ thesis/results/complete/imputed/pmm/results_pmm_", names[i], ".rds", sep = "")), ~.x$results)
}
################
# Loading data #
################
simdatasets <- list()
for (i in seq_len(nrow(combinations))) {
  simdatasets[[i]] <- map(read_rds(paste("/Volumes/Heleen\ 480GB/Master\ thesis/data/complete/simdata_", names[i], ".rds", sep = "")), ~.x %>% select(beta0j, beta1j, beta2j, beta3j, beta4j, beta5j, beta6j, beta7j, u0, u1, u2, u3, eij))
}
#####################################
# Creating functions for evaluation # 
#####################################
# Bias function
bias <- function(estimated, true) {
  estimates <- estimated %>%
    map(~ .x %>%
      filter(term == "(Intercept)" | term == "x1" | term == "x2" | term == "x3" | term == "x4" | term == "x5" | term == "x6" | term == "x7" | term == "z1" | term == "z2" | term == "x1:z1" | term == "x2:z1" | term == "x3:z2" | term == "sd__(Intercept)" | term == "sd__x1" | term == "sd__x2" | term == "sd__x3" | term == "sd__Observation") %>%
      pull(estimate, term))

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

  bias.datasets <- map2(estimates, truth, ~ .x - .y) %>% list_rbind() %>% as_tibble()
  bias <- colMeans(bias.datasets) %>% t() %>% as_tibble()

  return(list(bias.datasets = bias.datasets, bias = bias))
}
# MSE function
mse <- function(estimated, true) {
  estimates <- estimated %>%
    map(~ .x %>%
      filter(term == "(Intercept)" | term == "x1" | term == "x2" | term == "x3" | term == "x4" | term == "x5" | term == "x6" | term == "x7" | term == "z1" | term == "z2" | term == "x1:z1" | term == "x2:z1" | term == "x3:z2" | term == "sd__(Intercept)" | term == "sd__x1" | term == "sd__x2" | term == "sd__x3" | term == "sd__Observation") %>%
      pull(estimate, term))

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
bias.datasets_pmm <- map2(results_pmm, simdatasets, ~ bias(.x, .y))
bias_pmm <- map(bias.datasets_pmm, ~ .x$bias) %>% list_rbind()
# Saving
write_rds(bias.datasets_pmm, file = "/Volumes/Heleen\ 480GB/Master\ thesis/results/evaluations/bias_datasets_pmm.rds")
write_rds(bias_pmm, file = "/Volumes/Heleen\ 480GB/Master\ thesis/results/evaluations/bias_pmm.rds")
# MSE 
mse.datasets_pmm <- map2(results_pmm, simdatasets, ~ mse(.x, .y))
mse_pmm <- map(mse.datasets_pmm, ~ .x$mse) %>% list_rbind()
# Saving
write_rds(mse.datasets_pmm, file = "/Volumes/Heleen\ 480GB/Master\ thesis/results/evaluations/mse_datasets_pmm.rds")
write_rds(mse_pmm, file = "/Volumes/Heleen\ 480GB/Master\ thesis/results/evaluations/mse_pmm.rds")
# Coverage
coverage.datasets_pmm <- map2(results_pmm, simdatasets, ~ coverage(.x, .y))
coverage_pmm <- map(coverage.datasets_pmm, ~ .x$coverage) %>% list_rbind()
# Saving
write_rds(coverage.datasets_pmm, file = "/Volumes/Heleen\ 480GB/Master\ thesis/results/evaluations/coverage_datasets_pmm.rds")
write_rds(coverage_pmm, file = "/Volumes/Heleen\ 480GB/Master\ thesis/results/evaluations/coverage_pmm.rds")
#################
# Visualization #
#################
# Bias
bias_pmm <- read_rds("/Volumes/Heleen\ 480GB/Master\ thesis/results/evaluations/bias_pmm.rds")
bias_pmm <- cbind(bias_pmm, combinations) %>%
  pivot_longer(cols = beta0j:eij, names_to = "term", values_to = "bias") %>%
  mutate(ngroup = as.factor(ngroup), groupsize = as.factor(groupsize), icc = as.factor(icc), mar_mcar = as.factor(mar_mcar), miss = as.factor(miss), g = as.factor(g))

ggplot(bias_pmm, aes(
  x = bias,
  y = term,
  color = icc,
  shape = miss,
)) +
  geom_jitter() +
  facet_grid(rows = vars(ngroup, mar_mcar), cols = vars(groupsize, g), labeller = "label_both")  +
  geom_vline(xintercept = 0, color = "gray40") +
  geom_vline(xintercept = .10, linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = -.10, linetype = "dashed", color = "gray40") + 
  scale_color_viridis_d() + 
  theme_minimal() + 
  theme(panel.border = element_rect(colour = "gray25", fill = NA, size = .5), axis.text.x = element_text(size = 8, angle = 45, hjust = 1), legend.position = 'bottom') +
  labs(title = "Bias",
       x = "Bias",
       y = "Term",
       color = "ICC",
       shape = "Percentage of missing data")

# MSE 
mse_pmm <- read_rds("/Volumes/Heleen\ 480GB/Master\ thesis/results/evaluations/mse_pmm.rds")
mse_pmm <- cbind(mse_pmm, combinations) %>%
  pivot_longer(cols = beta0j:eij, names_to = "term", values_to = "mse") %>%
  mutate(ngroup = as.factor(ngroup), groupsize = as.factor(groupsize), icc = as.factor(icc), mar_mcar = as.factor(mar_mcar), miss = as.factor(miss), g = as.factor(g))

ggplot(mse_pmm, aes(
  x = mse,
  y = term,
  color = icc,
  shape = miss
)) +
  geom_jitter() +
  facet_grid(rows = vars(ngroup, mar_mcar), cols = vars(groupsize, g), labeller = "label_both")  + 
  scale_color_viridis_d() + 
  theme_minimal() + 
  theme(panel.border = element_rect(colour = "gray25", fill = NA, size = .5), axis.text.x = element_text(size = 8, angle = 45, hjust = 1), legend.position = 'bottom') +
  labs(title = "Bias",
       x = "Bias",
       y = "Term",
       color = "ICC",
       shape = "Percentage of missing data")

# Coverage
coverage_pmm <- read_rds("/Volumes/Heleen\ 480GB/Master\ thesis/results/evaluations/coverage_pmm.rds")
coverage_pmm <- cbind(coverage_pmm, combinations) %>%
  pivot_longer(cols = beta0j:`x3:z2`, names_to = "term", values_to = "coverage") %>%
  mutate(ngroup = as.factor(ngroup), groupsize = as.factor(groupsize), icc = as.factor(icc), mar_mcar = as.factor(mar_mcar), miss = as.factor(miss), g = as.factor(g))

ggplot(coverage_pmm, aes(
  x = coverage,
  y = term,
  color = icc,
  shape = miss
)) +
  geom_jitter() +
  facet_grid(rows = vars(ngroup, mar_mcar), cols = vars(groupsize, g), labeller = "label_both") +
  geom_vline(xintercept = .925, linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = .975, linetype = "dashed", color = "gray40") + 
  scale_color_viridis_d() + 
  theme_minimal() + 
  theme(panel.border = element_rect(colour = "gray25", fill = NA, size = .5), axis.text.x = element_text(size = 8, angle = 45, hjust = 1), legend.position = 'bottom') +
  labs(title = "Bias",
       x = "Bias",
       y = "Term",
       color = "ICC",
       shape = "Percentage of missing data")