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
# Source formatting #
#####################
source("scripts/Formatting for visualizations.R", encoding = "UTF-8")
#########
# Plots # 
#########
# Define color palette
cbbPalette <- c("#000000", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00", "#5D478B", "#0072B2", "#F0E442")
########
# Bias # 
########
# eij and u0
bias2 <- bias_combined %>%
    filter(Term == "eij" | Term == "u0") %>%
    ggplot(aes(
        x = Bias,
        y = Term,
        color = Method
    )) +
    # Absolute bias 
    geom_point(position = position_jitter(seed = 123, width = 0, height = .2)) + 
    # MCSE
    geom_errorbar(aes(xmin = Bias - MCSE, xmax = Bias + MCSE), position = position_jitter(seed = 123, width = 0, height = .2), width = .25) + 
    # Design factors
    facet_grid(cols = vars(`Number of groups`, `Group size`), rows = vars(`Missingness mechanism`), labeller = labeller(.rows = label_value, .cols = label_both)) +
    # Scale breaks for x-axis
    scale_x_continuous(n.breaks = 15, minor_breaks = seq(-40, 60, 2.5)) +
    # Scale labels y-axis for terms
    scale_y_discrete(limits = c("eij", "u0"), labels = c(expression(paste(epsilon, "ij")), expression(paste(upsilon, "0")))) +
    # Vertical line at 0
    geom_vline(xintercept = 0, color = "gray40") +
    # 10% relative bias lines
    geom_segment(y = 0, yend = 1.5, x = -2.5, linetype = "dashed", color = "gray40") +
    geom_segment(y = 0, yend = 1.5, x = 2.5, linetype = "dashed", color = "gray40") +
    geom_segment(y = 1.5, yend = 3, x = -8.474903, linetype = "dashed", color = "gray40") +
    geom_segment(y = 1.5, yend = 3, x = 8.474903, linetype = "dashed", color = "gray40") +
    # Color palette
    scale_color_manual(values = cbbPalette) +
    # Theme
    theme_minimal() +
    # Theme elements
    theme(panel.border = element_rect(colour = "gray25", fill = NA, size = .5), axis.text.x = element_text(size = 9, angle = 45, hjust = 1), axis.text.y = element_text(size = 9), legend.position = "bottom", legend.text = element_text(size = 9), axis.title = element_text(size = 12), legend.title = element_text(size = 12), panel.grid = element_line(color = "gray80")) +
    # Labels
    labs(
        x = "Absolute bias",
        y = "Term",
        color = "Method"
    )
# Saving plot 
ggsave(bias2, file = "thesis/graphs/bias2.png", width = 23, height = 9, units = "cm", dpi = "retina", bg = "white")
# Intercept 
biasintercept <- bias_combined %>%
    filter(Term == "beta0j") %>%
    ggplot(aes(
        x = Bias,
        y = Term,
        color = Method
    )) +
    # Absolute bias 
    geom_point(position = position_jitter(seed = 123, width = 0, height = .3)) +
    # MCSE 
    geom_errorbar(aes(xmin = Bias - MCSE, xmax = Bias + MCSE), position = position_jitter(seed = 123, width = 0, height = .3), width = .08) +
    # Design factors
    facet_grid(cols = vars(`Number of groups`, `Group size`), rows = vars(`Missingness mechanism`), labeller = labeller(.rows = label_value, .cols = label_both)) +
    # Scale breaks for x-axis
    scale_x_continuous(n.breaks = 10, minor_breaks = seq(-.9, .7, .05)) +
    # Scale labels y-axis for terms
    scale_y_discrete(limits = c("beta0j"), labels = c(expression(paste(gamma, "00")))) + 
    # Vertical line at 0
    geom_vline(xintercept = 0, color = "gray40") +
    # 10% relative bias lines
    geom_segment(y = 0, yend = 2, x = -.1, linetype = "dashed", color = "gray40") +
    geom_segment(y = 0, yend = 2, x = .1, linetype = "dashed", color = "gray40") +
    # Color palette
    scale_color_manual(values = cbbPalette) +
    # Theme
    theme_minimal() +
    # Theme elements
    theme(panel.border = element_rect(colour = "gray25", fill = NA, size = .5), axis.text.x = element_text(size = 9, angle = 45, hjust = 1), axis.text.y = element_text(size = 9), legend.position = "bottom", legend.text = element_text(size = 9), axis.title = element_text(size = 12), legend.title = element_text(size = 12), panel.grid = element_line(color = "gray80")) +
    labs(
        x = "Absolute bias",
        y = "Term",
        color = "Method"
    )
ggsave(biasintercept, file = "thesis/graphs/biasintercept.png",  width = 23, height = 9, units = "cm", dpi = "retina", bg = "white")
# Fixed level-1 effects 
biaslevel1 <- bias_combined %>%
    filter(Term == "beta1j" | Term == "beta2j" | Term == "beta3j" | Term == "beta4j" | Term == "beta5j" | Term == "beta6j" | Term == "beta7j") %>%
    ggplot(aes(
        x = Bias,
        y = Term,
        color = Method
    )) +
    # Absolute bias
    geom_point(position = position_jitter(seed = 123, width = 0, height = .2)) +
    # MCSE
    geom_errorbar(aes(xmin = Bias - MCSE, xmax = Bias + MCSE), position = position_jitter(seed = 123, width = 0, height = .2), width = .3) +
    # Design factors
    facet_grid(cols = vars(`Number of groups`, `Group size`), rows = vars(`Missingness mechanism`), labeller = labeller(.rows = label_value, .cols = label_both)) +
    # Scale breaks for x-axis
    scale_x_continuous(n.breaks = 10) +
    # Scale labels y-axis for terms
    scale_y_discrete(limits = c("beta1j", "beta2j", "beta3j", "beta4j", "beta5j", "beta6j", "beta7j"), labels = c(expression(paste(gamma, "10")), expression(paste(gamma, "20")), expression(paste(gamma, "30")), expression(paste(gamma, "40")), expression(paste(gamma, "50")), expression(paste(gamma, "60")), expression(paste(gamma, "70")))) + 
    # Vertical line at 0
    geom_vline(xintercept = 0, color = "gray40") +
    # 10% relative bias lines
    geom_vline(xintercept = -.05, linetype = "dashed", color = "gray40") +
    geom_vline(xintercept = .05, linetype = "dashed", color = "gray40") +
    # Color palette
    scale_color_manual(values = cbbPalette) +
    # Theme
    theme_minimal() +
    # Theme elements
    theme(panel.border = element_rect(colour = "gray25", fill = NA, size = .5), axis.text.x = element_text(size = 9, angle = 45, hjust = 1), axis.text.y = element_text(size = 9), legend.position = "bottom", legend.text = element_text(size = 9), axis.title = element_text(size = 12), legend.title = element_text(size = 12), panel.grid = element_line(color = "gray80")) +
    # Labels
    labs(
        x = "Absolute bias",
        y = "Term",
        color = "Method"
    )
# Saving plot
ggsave(biaslevel1, file = "thesis/graphs/biaslevel1.png",  width = 23, height = 15, units = "cm", dpi = "retina", bg = "white")
# Fixed level-2 effects
biaslevel2 <- bias_combined %>%
    filter(Term == "z1" | Term == "z2") %>%
    ggplot(aes(
        x = Bias,
        y = Term,
        color = Method
    )) +
    # Absolute bias
    geom_point(position = position_jitter(seed = 123, width = 0, height = .2)) +
    # MCSE
    geom_errorbar(aes(xmin = Bias - MCSE, xmax = Bias + MCSE), position = position_jitter(seed = 123, width = 0, height = .2), width = .25) +
    # Design factors
    facet_grid(cols = vars(`Number of groups`, `Group size`), rows = vars(`Missingness mechanism`), labeller = labeller(.rows = label_value, .cols = label_both)) +
    # Scale breaks for x-axis
    scale_x_continuous(n.breaks = 10, minor_breaks = seq(-.7, .9, .05)) +
    # Scale labels y-axis for terms
    scale_y_discrete(limits = c("z1", "z2"), labels = c(expression(paste(gamma, "01")), expression(paste(gamma, "02")))) + 
    # Vertical line at 0
    geom_vline(xintercept = 0, color = "gray40") +
    # 10% relative bias lines
    geom_segment(y = 0, yend = 3, x = -.05, linetype = "dashed", color = "gray40") +
    geom_segment(y = 0, yend = 3, x = .05, linetype = "dashed", color = "gray40") +
    # Color palette
    scale_color_manual(values = cbbPalette) +
    # Theme
    theme_minimal() +
    # Theme elements
    theme(panel.border = element_rect(colour = "gray25", fill = NA, size = .5), axis.text.x = element_text(size = 9, angle = 45, hjust = 1), axis.text.y = element_text(size = 9), legend.position = "bottom", legend.text = element_text(size = 9), axis.title = element_text(size = 12), legend.title = element_text(size = 12), panel.grid = element_line(color = "gray80")) +
    # Labels
    labs(
        x = "Absolute bias",
        y = "Term",
        color = "Method"
    )
# Saving plot
ggsave(biaslevel2, file = "thesis/graphs/biaslevel2.png",  width = 23, height = 9, units = "cm", dpi = "retina", bg = "white")
# Fixed cross-level effects
biascrosslevel <- bias_combined %>%
    filter(Term == "x1:z1" | Term == "x2:z1" | Term == "x3:z2") %>%
    ggplot(aes(
        x = Bias,
        y = Term,
        color = Method
    )) +
    # Absolute bias
    geom_point(position = position_jitter(seed = 123, width = 0, height = .1)) +
    # MCSE
    geom_errorbar(aes(xmin = Bias - MCSE, xmax = Bias + MCSE), position = position_jitter(seed = 123, width = 0, height = .1), width = .25) +
    # Design factors
    facet_grid(cols = vars(`Number of groups`, `Group size`), rows = vars(`Missingness mechanism`), labeller = labeller(.rows = label_value, .cols = label_both)) +
    # Scale breaks for x-axis
    scale_x_continuous(n.breaks = 10, minor_breaks = seq(-.325, .1, .025)) +
    # Scale labels y-axis for terms
    scale_y_discrete(limits = c("x1:z1", "x2:z1", "x3:z2"), labels = c(expression(paste(gamma, "11")), expression(paste(gamma, "21")), expression(paste(gamma, "32")))) + 
    # Vertical line at 0
    geom_vline(xintercept = 0, color = "gray40") +
    # 10% relative bias lines
    geom_segment(y = 0, yend = 4, x = -.035, linetype = "dashed", color = "gray40") +
    geom_segment(y = 0, yend = 4, x = .035, linetype = "dashed", color = "gray40") +
    # Color palette
    scale_color_manual(values = cbbPalette) +
    # Theme
    theme_minimal() +
    # Theme elements
    theme(panel.border = element_rect(colour = "gray25", fill = NA, size = .5), axis.text.x = element_text(size = 9, angle = 45, hjust = 1), axis.text.y = element_text(size = 9), legend.position = "bottom", legend.text = element_text(size = 9), axis.title = element_text(size = 12), legend.title = element_text(size = 12), panel.grid = element_line(color = "gray80")) +
    # Labels
    labs(
        x = "Absolute bias",
        y = "Term",
        color = "Method"
    )
# Saving plot
ggsave(biascrosslevel, file = "thesis/graphs/biascrosslevel.png",  width = 23, height = 10, units = "cm", dpi = "retina", bg = "white")
# Random slopes
biasrandom <- bias_combined %>%
    filter(Term == "u1" | Term == "u2" | Term == "u3") %>%
    ggplot(aes(
        x = Bias,
        y = Term,
        color = Method
    )) +
    # Absolute bias
    geom_point(position = position_jitter(seed = 123, width = 0, height = .1)) +
    # MCSE
    geom_errorbar(aes(xmin = Bias - MCSE, xmax = Bias + MCSE), position = position_jitter(seed = 123, width = 0, height = .1), width = .25) +
    # Design factors
    facet_grid(cols = vars(`Number of groups`, `Group size`), rows = vars(`Missingness mechanism`), labeller = labeller(.rows = label_value, .cols = label_both)) +
    # Scale breaks for x-axis
    scale_x_continuous(n.breaks = 10, minor_breaks = seq(-.6, .25, .025)) +
    # Scale labels y-axis for terms
    scale_y_discrete(limits = c("u1", "u2", "u3"), labels = c(expression(paste(upsilon, "1")), expression(paste(upsilon, "2")), expression(paste(upsilon, "3")))) +
    # Vertical line at 0 
    geom_vline(xintercept = 0, color = "gray40") +
    # 10% relative bias lines
    geom_segment(y = 0, yend = 4, x = -.1, linetype = "dashed", color = "gray40") +
    geom_segment(y = 0, yend = 4, x = .1, linetype = "dashed", color = "gray40") +
    # Color palette
    scale_color_manual(values = cbbPalette) +
    # Theme
    theme_minimal() +
    # Theme elements
    theme(panel.border = element_rect(colour = "gray25", fill = NA, size = .5), axis.text.x = element_text(size = 9, angle = 45, hjust = 1), axis.text.y = element_text(size = 9), legend.position = "bottom", legend.text = element_text(size = 9), axis.title = element_text(size = 12), legend.title = element_text(size = 12), panel.grid = element_line(color = "gray80")) +
    # Labels
    labs(
        x = "Absolute bias",
        y = "Term",
        color = "Method"
    )
# Saving plot
ggsave(biasrandom, file = "thesis/graphs/biasrandom.png",  width = 23, height = 10, units = "cm", dpi = "retina", bg = "white")
############
# Coverage #
############
# Fixed level-1 effects
coveragelevel1 <- coverage_combined %>%
    filter(Term == "beta0j" | Term == "beta1j" | Term == "beta2j" | Term == "beta3j" | Term == "beta4j" | Term == "beta5j" | Term == "beta6j" | Term == "beta7j") %>%
    ggplot(aes(
        x = Coverage,
        y = Term,
        color = Method
    )) +
    # Coverage
    geom_point(position = position_jitter(seed = 123, width = 0, height = .2)) +
    # MCSE
    geom_errorbar(aes(xmin = Coverage - MCSE, xmax = Coverage + MCSE), position = position_jitter(seed = 123, width = 0, height = .2), width = .3) +
    # Vertical line at 95% coverage
    geom_vline(xintercept = .95, color = "gray40") +
    # Dashed lines at 92.5% and 97.5% coverage
    geom_vline(xintercept = .925, linetype = "dashed", color = "gray40") +
    geom_vline(xintercept = .975, linetype = "dashed", color = "gray40") +
    # Design factors
    facet_grid(cols = vars(`Number of groups`, `Group size`), rows = vars(`Missingness mechanism`), labeller = labeller(.rows = label_value, .cols = label_both)) +
    # Scale breaks for x-axis
    scale_x_continuous(n.breaks = 10, minor_breaks = seq(0, 1, .025)) +
    # Scale labels y-axis for terms
    scale_y_discrete(limits = c("beta0j", "beta1j", "beta2j", "beta3j", "beta4j", "beta5j", "beta6j", "beta7j"), labels = c(expression(paste(gamma, "00")), expression(paste(gamma, "10")), expression(paste(gamma, "20")), expression(paste(gamma, "30")), expression(paste(gamma, "40")), expression(paste(gamma, "50")), expression(paste(gamma, "60")), expression(paste(gamma, "70")))) +
    # Color palette
    scale_color_manual(values = cbbPalette) +
    # Theme
    theme_minimal() +
    # Theme elements
    theme(panel.border = element_rect(colour = "gray25", fill = NA, size = .5), axis.text.x = element_text(size = 9, angle = 45, hjust = 1), axis.text.y = element_text(size = 9), legend.position = "bottom", legend.text = element_text(size = 9), axis.title = element_text(size = 12), legend.title = element_text(size = 12), panel.grid = element_line(color = "gray80")) +
    # Labels
    labs(
        x = "Coverage of 95% confidence intervals",
        y = "Term",
        color = "Method"
    )
# Saving plot
ggsave(coveragelevel1, file = "thesis/graphs/coveragelevel1.png", width = 23, height = 15, units = "cm", dpi = "retina", bg = "white")
# Fixed level-2 effects
coveragelevel2 <- coverage_combined %>%
    filter(Term == "z1" | Term == "z2" | Term == "x1:z1" | Term == "x2:z1" | Term == "x3:z2") %>%
    ggplot(aes(
        x = Coverage,
        y = Term,
        color = Method
    )) +
    # Coverage
    geom_point(position = position_jitter(seed = 123, width = 0, height = .2)) +
    # MCSE
    geom_errorbar(aes(xmin = Coverage - MCSE, xmax = Coverage + MCSE), position = position_jitter(seed = 123, width = 0, height = .2), width = .25) +
    # Vertical line at 95% coverage
    geom_vline(xintercept = .95, color = "gray40") +
    # Dashed lines at 92.5% and 97.5% coverage
    geom_vline(xintercept = .925, linetype = "dashed", color = "gray40") +
    geom_vline(xintercept = .975, linetype = "dashed", color = "gray40") +
    # Design factors
    facet_grid(cols = vars(`Number of groups`, `Group size`), rows = vars(`Missingness mechanism`), labeller = labeller(.rows = label_value, .cols = label_both)) +
    # Scale breaks for x-axis
    scale_x_continuous(n.breaks = 10, minor_breaks = seq(0, 1, .025)) +
    # Scale labels y-axis for terms
    scale_y_discrete(limits = c("z1", "z2", "x1:z1", "x2:z1", "x3:z2"), labels = c(expression(paste(gamma, "01")), expression(paste(gamma, "02")), expression(paste(gamma, "11")), expression(paste(gamma, "21")), expression(paste(gamma, "32")))) +
    # Color palette
    scale_color_manual(values = cbbPalette) +
    # Theme
    theme_minimal() +
    # Theme elements
    theme(panel.border = element_rect(colour = "gray25", fill = NA, size = .5), axis.text.x = element_text(size = 9, angle = 45, hjust = 1), axis.text.y = element_text(size = 9), legend.position = "bottom", legend.text = element_text(size = 9), axis.title = element_text(size = 12), legend.title = element_text(size = 12), panel.grid = element_line(color = "gray80")) +
    # Labels
    labs(
        x = "Coverage of 95% confidence intervals",
        y = "Term",
        color = "Method"
    )
# Saving plot
ggsave(coveragelevel2, file = "thesis/graphs/coveragelevel2.png", width = 23, height = 13, units = "cm", dpi = "retina", bg = "white")
#######
# CIW #
#######
# Intercept
ciwintercept <- ciw_combined %>%
    filter(Term == "beta0j") %>%
    ggplot(aes(
        x = CIW,
        y = Term,
        color = Method
    )) +
    # CIW
    geom_point(position = position_jitter(seed = 123, width = 0, height = .1)) +
    # Design factors
    facet_grid(cols = vars(`Number of groups`, `Group size`), rows = vars(`Missingness mechanism`), labeller = labeller(.rows = label_value, .cols = label_both)) +
    # Scale breaks for x-axis
    scale_x_continuous(n.breaks = 5, minor_breaks = seq(4, 8, .25)) +
    # Scale labels y-axis for terms
    scale_y_discrete(limits = c("beta0j"), labels = c(expression(paste(gamma, "00")))) +
    # Color palette
    scale_color_manual(values = cbbPalette) +
    # Theme
    theme_minimal() +
    # Theme elements
    theme(panel.border = element_rect(colour = "gray25", fill = NA, size = .5), axis.text.x = element_text(size = 9, angle = 0, hjust = 1), axis.text.y = element_text(size = 9), legend.position = "bottom", legend.text = element_text(size = 9), axis.title = element_text(size = 12), legend.title = element_text(size = 12), panel.grid = element_line(color = "gray80")) +
    # Labels
    labs(
        x = "Confidence interval width",
        y = "Term",
        color = "Method"
    )
# Saving plot
ggsave(ciwintercept, file = "thesis/graphs/ciwintercept.png", width = 23, height = 8, units = "cm", dpi = "retina", bg = "white")
# Fixed level-1 effects
ciwlevel1 <- ciw_combined %>%
    filter(Term == "beta1j" | Term == "beta2j" | Term == "beta3j" | Term == "beta4j" | Term == "beta5j" | Term == "beta6j" | Term == "beta7j") %>%
    ggplot(aes(
        x = CIW,
        y = Term,
        color = Method
    )) +
    # CIW
    geom_point(position = position_jitter(seed = 123, width = 0, height = .1)) +
    # Design factors
    facet_grid(cols = vars(`Number of groups`, `Group size`), rows = vars(`Missingness mechanism`), labeller = "label_both") +
    # Scale breaks for x-axis
    scale_x_continuous(n.breaks = 10) +
    # Scale labels y-axis for terms
    scale_y_discrete(limits = c("beta1j", "beta2j", "beta3j", "beta4j", "beta5j", "beta6j", "beta7j"), labels = c(expression(paste(gamma, "10")), expression(paste(gamma, "20")), expression(paste(gamma, "30")), expression(paste(gamma, "40")), expression(paste(gamma, "50")), expression(paste(gamma, "60")), expression(paste(gamma, "70")))) +
    # Color palette
    scale_color_manual(values = cbbPalette) +
    # Theme
    theme_minimal() +
    # Theme elements
    theme(panel.border = element_rect(colour = "gray25", fill = NA, size = .5), axis.text.x = element_text(size = 9, angle = 45, hjust = 1), axis.text.y = element_text(size = 9), legend.position = "bottom", legend.text = element_text(size = 9), axis.title = element_text(size = 12), legend.title = element_text(size = 12), panel.grid = element_line(color = "gray80")) +
    # Labels
    labs(
        x = "Confidence interval width",
        y = "Term",
        color = "Method"
    )
# Saving plot
ggsave(ciwlevel1, file = "thesis/graphs/ciwlevel1.png", width = 23, height = 15, units = "cm", dpi = "retina", bg = "white")
# Fixed level-2 effects
ciwlevel2 <- ciw_combined %>%
    filter(Term == "z1" | Term == "z2") %>%
    ggplot(aes(
        x = CIW,
        y = Term,
        color = Method
    )) +
    # CIW
    geom_point(position = position_jitter(seed = 123, width = 0, height = .2)) +
    # Design factors
    facet_grid(cols = vars(`Number of groups`, `Group size`), rows = vars(`Missingness mechanism`), labeller = labeller(.rows = label_value, .cols = label_both)) +
    # Scale breaks for x-axis
    scale_x_continuous(n.breaks = 10, minor_breaks = seq(0, 8, .25)) +
    # Scale labels y-axis for terms
    scale_y_discrete(limits = c("z1", "z2"), labels = c(expression(paste(gamma, "01")), expression(paste(gamma, "02")))) +
    # Color palette
    scale_color_manual(values = cbbPalette) +
    # Theme
    theme_minimal() +
    # Theme elements
    theme(panel.border = element_rect(colour = "gray25", fill = NA, size = .5), axis.text.x = element_text(size = 9, angle = 0, hjust = 1), axis.text.y = element_text(size = 9), legend.position = "bottom", legend.text = element_text(size = 9), axis.title = element_text(size = 12), legend.title = element_text(size = 12), panel.grid = element_line(color = "gray80")) +
    # Labels
    labs(
        x = "Confidence interval width",
        y = "Term",
        color = "Method"
    )
# Saving plot
ggsave(ciwlevel2, file = "thesis/graphs/ciwlevel2.png", width = 23, height = 9, units = "cm", dpi = "retina", bg = "white")
# Cross-level effects
ciwcrosslevel <- ciw_combined %>%
    filter(Term == "x1:z1" | Term == "x2:z1" | Term == "x3:z2") %>%
    ggplot(aes(
        x = CIW,
        y = Term,
        color = Method
    )) +
    # CIW
    geom_point(position = position_jitter(seed = 123, width = 0, height = .2)) +
    # Design factors
    facet_grid(cols = vars(`Number of groups`, `Group size`), rows = vars(`Missingness mechanism`), labeller = labeller(.rows = label_value, .cols = label_both)) +
    # Scale breaks for x-axis
    scale_x_continuous(n.breaks = 10) +
    # Scale labels y-axis for terms
    scale_y_discrete(limits = c("x1:z1", "x2:z1", "x3:z2"), labels = c(expression(paste(gamma, "11")), expression(paste(gamma, "21")), expression(paste(gamma, "32")))) +
    # Color palette
    scale_color_manual(values = cbbPalette) +
    # Theme
    theme_minimal() +
    # Theme elements
    theme(panel.border = element_rect(colour = "gray25", fill = NA, size = .5), axis.text.x = element_text(size = 9, angle = 45, hjust = 1), axis.text.y = element_text(size = 9), legend.position = "bottom", legend.text = element_text(size = 9), axis.title = element_text(size = 12), legend.title = element_text(size = 12), panel.grid = element_line(color = "gray80")) +
    # Labels
    labs(
        x = "Confidence interval width",
        y = "Term",
        color = "Method"
    )
# Saving plot
ggsave(ciwcrosslevel, file = "thesis/graphs/ciwcrosslevel.png", width = 23, height = 10, units = "cm", dpi = "retina", bg = "white")
