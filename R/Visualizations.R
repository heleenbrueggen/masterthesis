##################
# Visualizations #
##################
#############
# Libraries #
#############
library(ggplot2)
library(tidyverse)
################
# Setting seed #
################
set.seed(123)
#####################
# Complete analysis #
#####################
#######################
# Defining parameters #
#######################
ngroups <- c(30, 50)
groupsizes <- c(15, 35, 50)
iccs <- c(.2, .5)
mar_mcar <- c("mcar")
miss <- c(0)
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
# Bias
bias_nomiss <- read_rds("/Volumes/Heleen\ 480GB/Master\ thesis/results/evaluations/bias_nomiss.rds")
bias_nomiss <- cbind(bias_nomiss, combinations) %>%
  pivot_longer(cols = beta0j:eij, names_to = "term", values_to = "bias") %>%
  mutate(ngroup = as.factor(ngroup), groupsize = as.factor(groupsize), icc = as.factor(icc), mar_mcar = as.factor(mar_mcar), miss = as.factor(miss), g = as.factor(g))

# MSE 
mse_nomiss <- read_rds("/Volumes/Heleen\ 480GB/Master\ thesis/results/evaluations/mse_nomiss.rds")
mse_nomiss <- cbind(mse_nomiss, combinations) %>%
  pivot_longer(cols = beta0j:eij, names_to = "term", values_to = "mse") %>%
  mutate(ngroup = as.factor(ngroup), groupsize = as.factor(groupsize), icc = as.factor(icc), mar_mcar = as.factor(mar_mcar), miss = as.factor(miss), g = as.factor(g))

# Coverage
coverage_nomiss <- read_rds("/Volumes/Heleen\ 480GB/Master\ thesis/results/evaluations/coverage_nomiss.rds")
coverage_nomiss <- cbind(coverage_nomiss, combinations) %>%
  pivot_longer(cols = beta0j:`x3:z2`, names_to = "term", values_to = "coverage") %>%
  mutate(ngroup = as.factor(ngroup), groupsize = as.factor(groupsize), icc = as.factor(icc), mar_mcar = as.factor(mar_mcar), miss = as.factor(miss), g = as.factor(g))
#######
# PMM # 
#######
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
# Bias 
bias_pmm <- read_rds("/Volumes/Heleen\ 480GB/Master\ thesis/results/evaluations/bias_pmm.rds")
bias_pmm <- cbind(bias_pmm, combinations) %>%
  pivot_longer(cols = beta0j:eij, names_to = "term", values_to = "bias") %>%
  mutate(ngroup = as.factor(ngroup), groupsize = as.factor(groupsize), icc = as.factor(icc), mar_mcar = as.factor(mar_mcar), miss = as.factor(miss), g = as.factor(g))

# MSE 
mse_pmm <- read_rds("/Volumes/Heleen\ 480GB/Master\ thesis/results/evaluations/mse_pmm.rds")
mse_pmm <- cbind(mse_pmm, combinations) %>%
  pivot_longer(cols = beta0j:eij, names_to = "term", values_to = "mse") %>%
  mutate(ngroup = as.factor(ngroup), groupsize = as.factor(groupsize), icc = as.factor(icc), mar_mcar = as.factor(mar_mcar), miss = as.factor(miss), g = as.factor(g))

# Coverage
coverage_pmm <- read_rds("/Volumes/Heleen\ 480GB/Master\ thesis/results/evaluations/coverage_pmm.rds")
coverage_pmm <- cbind(coverage_pmm, combinations) %>%
  pivot_longer(cols = beta0j:`x3:z2`, names_to = "term", values_to = "coverage") %>%
  mutate(ngroup = as.factor(ngroup), groupsize = as.factor(groupsize), icc = as.factor(icc), mar_mcar = as.factor(mar_mcar), miss = as.factor(miss), g = as.factor(g))
#####################
# Combining results # 
#####################
bias_combined <- rbind(
    bias_nomiss %>% mutate(method = "complete"), 
    bias_pmm %>% mutate(method = "pmm"))
colnames(bias_combined) <- c("Number of groups", "Group size", "ICC", "Missingness mechanism", "Percentage of missing data", "Gamma", "Term", "Bias", "Method")
mse_combined <- rbind(
    mse_nomiss %>% mutate(method = "complete"), 
    mse_pmm %>% mutate(method = "pmm"))
colnames(mse_combined) <- c("Number of groups", "Group size", "ICC", "Missingness mechanism", "Percentage of missing data", "Gamma", "Term", "MSE", "Method")
coverage_combined <- rbind(
    coverage_nomiss %>% mutate(method = "complete"), 
    coverage_pmm %>% mutate(method = "pmm"))
colnames(coverage_combined) <- c("Number of groups", "Group size", "ICC", "Missingness mechanism", "Percentage of missing data", "Gamma", "Term", "Coverage", "Method")
#########
# Plots # 
#########
# Bias 
bias_combined %>%
    filter(Method == "complete") %>%
    ggplot(aes(
        x = Bias,
        y = Term,
        color = ICC
    )) +
    geom_jitter(size = 2) +
    facet_grid(rows = vars(Gamma), cols = vars(`Group size`, `Number of groups`), labeller = "label_both", scales = "free_x") +
    geom_vline(xintercept = 0, color = "gray40") +
    geom_vline(xintercept = .10, linetype = "dashed", color = "gray40") +
    geom_vline(xintercept = -.10, linetype = "dashed", color = "gray40") +
    # scale_color_viridis_d() +
    scale_color_manual(values = c("firebrick", "chartreuse3")) +
    theme_minimal() +
    theme(panel.border = element_rect(colour = "gray25", fill = NA, size = .5), axis.text.x = element_text(size = 8, angle = 45, hjust = 1), legend.position = "bottom") +
    labs(
        title = "Bias for Complete",
        x = "Bias",
        y = "Term",
        color = "ICC"
    )
bias_combined %>%
    filter(Method == "pmm") %>%
    ggplot(aes(
        x = Bias,
        y = Term,
        color = ICC,
        shape = `Percentage of missing data`
    )) +
    geom_jitter(size = 2) +
    facet_grid(rows = vars(Gamma, `Missingness mechanism`), cols = vars(`Group size`, `Number of groups`), labeller = "label_both", scales = "free_x") +
    geom_vline(xintercept = 0, color = "gray40") +
    geom_vline(xintercept = .10, linetype = "dashed", color = "gray40") +
    geom_vline(xintercept = -.10, linetype = "dashed", color = "gray40") +
    scale_shape_manual(values = c(3, 20)) +
    # scale_color_viridis_d() +
    scale_color_manual(values = c("firebrick", "chartreuse3")) +
    theme_minimal() +
    theme(panel.border = element_rect(colour = "gray25", fill = NA, size = .5), axis.text.x = element_text(size = 8, angle = 45, hjust = 1), legend.position = "bottom") +
    labs(
        title = "Bias for PMM",
        x = "Bias",
        y = "Term",
        color = "ICC",
        shape = "Percentage of missing data"
    )

# MSE 
mse_combined %>%
    filter(Method == "complete") %>%
    ggplot(aes(
        x = MSE,
        y = Term,
        color = ICC
    )) +
    geom_jitter() +
    facet_grid(rows = vars(Gamma), cols = vars(`Group size`, `Number of groups`), labeller = "label_both") +
    # scale_color_viridis_d() +
    scale_color_manual(values = c("firebrick", "chartreuse3")) +
    theme_minimal() +
    theme(panel.border = element_rect(colour = "gray25", fill = NA, size = .5), axis.text.x = element_text(size = 8, hjust = 1), legend.position = "bottom") +
    labs(
        title = "MSE for Complete",
        x = "MSE",
        y = "Term",
        color = "ICC"
    )
mse_combined %>%
    filter(Method == "pmm") %>%
    ggplot(aes(
        x = MSE,
        y = Term,
        color = ICC,
        shape = `Percentage of missing data`
    )) +
    geom_jitter() +
    facet_grid(rows = vars(Gamma, `Missingness mechanism`), cols = vars(`Group size`, `Number of groups`), labeller = "label_both") +
    scale_shape_manual(values = c(3, 20)) +
    # scale_color_viridis_d() +
    scale_color_manual(values = c("firebrick", "chartreuse3")) +
    theme_minimal() +
    theme(panel.border = element_rect(colour = "gray25", fill = NA, size = .5), axis.text.x = element_text(size = 8, angle = 45, hjust = 1), legend.position = "bottom") +
    labs(
        title = "MSE for PMM",
        x = "MSE",
        y = "Term",
        color = "ICC",
        shape = "Percentage of missing data"
    )
# Coverage
coverage_combined %>%
    filter(Method == "complete") %>%
    ggplot(aes(
        x = Coverage,
        y = Term,
        color = ICC
    )) +
    geom_jitter() +
    facet_grid(rows = vars(Gamma), cols = vars(`Group size`, `Number of groups`), labeller = "label_both") +
    geom_vline(xintercept = .925, linetype = "dashed", color = "gray40") +
    geom_vline(xintercept = .975, linetype = "dashed", color = "gray40") +
    #   scale_color_viridis_d() +
    scale_color_manual(values = c("firebrick", "chartreuse3")) +
    theme_minimal() +
    theme(panel.border = element_rect(colour = "gray25", fill = NA, size = .5), axis.text.x = element_text(size = 8, angle = 45, hjust = 1), legend.position = "bottom") +
    labs(
        title = "Coverage for Complete",
        x = "Coverage",
        y = "Term",
        color = "ICC"
    )
coverage_combined %>%
    filter(Method == "pmm") %>%
    ggplot(aes(
        x = Coverage,
        y = Term,
        color = ICC,
        shape = `Percentage of missing data`
    )) +
    geom_jitter() +
    facet_grid(rows = vars(Gamma, `Missingness mechanism`), cols = vars(`Group size`, `Number of groups`), labeller = "label_both") +
    geom_vline(xintercept = .925, linetype = "dashed", color = "gray40") +
    geom_vline(xintercept = .975, linetype = "dashed", color = "gray40") +
    scale_shape_manual(values = c(3, 20)) +
    #   scale_color_viridis_d() +
    scale_color_manual(values = c("firebrick", "chartreuse3")) +
    theme_minimal() +
    theme(panel.border = element_rect(colour = "gray25", fill = NA, size = .5), axis.text.x = element_text(size = 8, angle = 45, hjust = 1), legend.position = "bottom") +
    labs(
        title = "Coverage for PMM",
        x = "Coverage",
        y = "Term",
        color = "ICC",
        shape = "Percentage of missing data"
    )










bias_combined %>%
    ggplot(aes(
        x = bias,
        y = term,
        color = icc,
        shape = miss
    )) +
    geom_jitter(size = 2) +
    facet_grid(rows = vars(ngroup, groupsize), cols = vars(method, mar_mcar), labeller = "label_both", scales = "free_x") +
    geom_vline(xintercept = 0, color = "gray40") +
    geom_vline(xintercept = .10, linetype = "dashed", color = "gray40") +
    geom_vline(xintercept = -.10, linetype = "dashed", color = "gray40") +
    scale_shape_manual(values = c(0, 3, 20)) +
    # scale_color_viridis_d() +
    scale_color_manual(values = c("firebrick", "chartreuse3")) +
    theme_minimal() +
    theme(panel.border = element_rect(colour = "gray25", fill = NA, size = .5), axis.text.x = element_text(size = 8, angle = 45, hjust = 1), legend.position = "bottom") +
    labs(
        title = "Bias",
        x = "Bias",
        y = "Term",
        color = "ICC",
        shape = "Percentage of missing data"
    )
# beta0j
bias_combined %>%
    filter(term == "beta0j") %>%
    ggplot(aes(
        x = bias,
        y = method,
        color = icc,
        shape = miss
    )) +
    geom_jitter(size = 2) +
    facet_grid(rows = vars(ngroup, groupsize), cols = vars(mar_mcar, g), labeller = "label_both", scales = "free_x") +
    geom_vline(xintercept = 0, color = "gray40") +
    geom_vline(xintercept = .10, linetype = "dashed", color = "gray40") +
    geom_vline(xintercept = -.10, linetype = "dashed", color = "gray40") +
    scale_shape_manual(values = c(0, 3, 20)) +
    # scale_color_viridis_d() +
    scale_color_manual(values = c("firebrick", "chartreuse3")) +
    theme_minimal() +
    theme(panel.border = element_rect(colour = "gray25", fill = NA, size = .5), axis.text.x = element_text(size = 8, angle = 45, hjust = 1), legend.position = "bottom") +
    labs(
        title = "Bias of beta0j",
        x = "Bias",
        y = "Term",
        color = "ICC",
        shape = "Percentage of missing data"
    )