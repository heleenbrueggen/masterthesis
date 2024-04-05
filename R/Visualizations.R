##################
# Visualizations #
##################
#############
# Libraries #
#############
library(ggplot2)
library(ggsci)
library(tidyverse)
################
# Setting seed #
################
set.seed(123)
#####################
# Source formatting #
#####################
source("/Users/Heleen/Desktop/Universiteit Utrecht/Methodology & Statistics for the Behavioural, Biomedical and Social Sciences/Year 2/Master Thesis/masterthesis/R/Formatting for visualizations.R", encoding = "UTF-8")
#########
# Plots # 
#########
# Bias 
bias_combined %>%
    filter(Term != "eij" & Term != "u0") %>%
    ggplot(aes(
        x = Bias,
        y = Term,
        color = Method
    )) +
    geom_point(position = position_jitter(seed = 123), size = 2) +
    geom_errorbar(aes(xmin = Bias - MCSE, xmax = Bias + MCSE), position = position_jitter(seed = 123), width = .2) +
    # geom_pointrange(aes(xmin = Bias - MCSE, xmax = Bias + MCSE), position = position_jitter(seed = 123)) + 
    # geom_linerange(aes(xmin = Bias - MCSE, xmax = Bias + MCSE), position = position_jitter(seed = 123)) +  
    facet_grid(cols = vars(`Number of groups`, `Group size`), rows = vars(`Missingness mechanism`), labeller = "label_both") +
    scale_x_continuous(limits = c(-1, 1), breaks = c(-1, -.75, -.5, -.25, 0, .25, .5, .75, 1)) +
    geom_vline(xintercept = 0, color = "gray40") +
    # geom_vline(xintercept = .10, linetype = "dashed", color = "gray40") +
    # geom_vline(xintercept = -.10, linetype = "dashed", color = "gray40") +
    # scale_color_viridis_d() +
    # scale_color_brewer(palette = "Dark2") +
    # scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#E6AB02", "#666666")) +
    scale_color_lancet() +
    theme_minimal() +
    theme(panel.border = element_rect(colour = "gray25", fill = NA, size = .5), axis.text.x = element_text(size = 8, angle = 45, hjust = 1), legend.position = "bottom") +
    labs(
        title = "Bias",
        x = "Bias",
        y = "Term",
        color = "Method"
    )
bias_combined %>%
    filter(Term == "eij" | Term == "u0") %>%
    ggplot(aes(
        x = Bias,
        y = Term,
        color = Method
    )) +
    geom_point(position = position_jitter(seed = 123), size = 2) +
    geom_errorbar(aes(xmin = Bias - MCSE, xmax = Bias + MCSE), position = position_jitter(seed = 123), width = .02) +
    facet_grid(cols = vars(`Number of groups`, `Group size`), rows = vars(`Missingness mechanism`), labeller = "label_both") +
    geom_vline(xintercept = 0, color = "gray40") +
    # geom_vline(xintercept = .10, linetype = "dashed", color = "gray40") +
    # geom_vline(xintercept = -.10, linetype = "dashed", color = "gray40") +
    scale_color_lancet() +
    theme_minimal() +
    theme(panel.border = element_rect(colour = "gray25", fill = NA, size = .5), axis.text.x = element_text(size = 8, angle = 45, hjust = 1), legend.position = "bottom") +
    labs(
        title = "Bias",
        x = "Bias",
        y = "Term",
        color = "Method"
    )
# MSE
mse_combined %>%
    filter(Term != "eij" & Term != "u0") %>%
    ggplot(aes(
        x = MSE,
        y = Term,
        color = Method
    )) +
    geom_jitter() +
    facet_grid(cols = vars(`Number of groups`, `Group size`), rows = vars(`Missingness mechanism`), labeller = "label_both") +
    scale_color_lancet() +
    theme_minimal() +
    theme(panel.border = element_rect(colour = "gray25", fill = NA, size = .5), axis.text.x = element_text(size = 8, hjust = 1), legend.position = "bottom") +
    labs(
        title = "MSE",
        x = "MSE",
        y = "Term",
        color = "Method"
    ) 
mse_combined %>%
    filter(Term == "eij" | Term == "u0") %>%
    ggplot(aes(
        x = MSE,
        y = Term,
        color = Method
    )) +
    geom_jitter() +
    facet_grid(cols = vars(`Number of groups`, `Group size`), rows = vars(`Missingness mechanism`), labeller = "label_both") +
    scale_color_lancet() +
    theme_minimal() +
    theme(panel.border = element_rect(colour = "gray25", fill = NA, size = .5), axis.text.x = element_text(size = 8, hjust = 1), legend.position = "bottom") +
    labs(
        title = "MSE",
        x = "MSE",
        y = "Term",
        color = "Method"
    ) 
# Coverage
coverage_combined %>%
    ggplot(aes(
        x = Coverage,
        y = Term,
        color = Method
    )) +
    geom_point(position = position_jitter(seed = 123), size = 2) +
    geom_errorbar(aes(xmin = Coverage - MCSE, xmax = Coverage + MCSE), position = position_jitter(seed = 123), width = .2) +
    geom_vline(xintercept = .925, linetype = "dashed", color = "gray40") +
    geom_vline(xintercept = .975, linetype = "dashed", color = "gray40") +
    facet_grid(cols = vars(`Number of groups`, `Group size`), rows = vars(`Missingness mechanism`), labeller = "label_both") +
    scale_color_lancet() +
    theme_minimal() +
    theme(panel.border = element_rect(colour = "gray25", fill = NA, size = .5), axis.text.x = element_text(size = 8, angle = 45, hjust = 1), legend.position = "bottom") +
    labs(
        title = "Coverage",
        x = "Coverage",
        y = "Term",
        color = "Method"
    )
# CIW 
ciw_combined %>%
    ggplot(aes(
        x = CIW,
        y = Term,
        color = Method
    )) +
    geom_point(position = position_jitter(seed = 123), size = 2) +
    # geom_errorbar(aes(xmin = Coverage - MCSE, xmax = Coverage + MCSE), position = position_jitter(seed = 123), width = .2) +
    # geom_linerange(aes(xmin = 0, xmax = CIW), position = position_jitter(seed = 123)) +
    facet_grid(cols = vars(`Number of groups`, `Group size`), rows = vars(`Missingness mechanism`), labeller = "label_both") +
    scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8)) +
    scale_color_lancet() +
    theme_minimal() +
    theme(panel.border = element_rect(colour = "gray25", fill = NA, size = .5), axis.text.x = element_text(size = 8, hjust = 1), legend.position = "bottom") +
    labs(
        title = "Confidence interval width",
        x = "CIW",
        y = "Term",
        color = "Method"
    )
