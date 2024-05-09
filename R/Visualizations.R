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
cbbPalette <- c("#000000", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00", "#5D478B", "#0072B2", "#F0E442")
# cPalette <- c("#e60049", "#0bb4ff", "#50e991", "#e6d800", "#9b19f5", "#ffa300", "#dc0ab4", "#b3d4ff", "#00bfa0")
# cPalette <- c("#fd7f6f", "#7eb0d5", "#b2e061", "#bd7ebe", "#ffb55a", "#ffee65", "#beb9db", "#fdcce5", "#8bd3c7")
########
# Bias # 
########
bias1 <- bias_combined %>%
    filter(Term != "eij" & Term != "u0") %>%
    ggplot(aes(
        x = Bias,
        y = Term,
        color = Method
    )) +
    geom_point(position = position_jitter(seed = 123, width = 0)) +
    geom_errorbar(aes(xmin = Bias - MCSE, xmax = Bias + MCSE), position = position_jitter(seed = 123, width = 0), width = .3) +
    # geom_pointrange(aes(xmin = Bias - MCSE, xmax = Bias + MCSE), position = position_jitter(seed = 123)) +
    # geom_linerange(aes(xmin = Bias - MCSE, xmax = Bias + MCSE), position = position_jitter(seed = 123)) +
    facet_grid(cols = vars(`Number of groups`, `Group size`), rows = vars(`Missingness mechanism`), labeller = labeller(.rows = label_value, .cols = label_both)) +
    scale_x_continuous(n.breaks = 20) +
    scale_y_discrete(limits = c("beta0j", "beta1j", "beta2j", "beta3j", "beta4j", "beta5j", "beta6j", "beta7j", "z1", "z2", "x1:z1", "x2:z1", "x3:z2", "u1", "u2", "u3"), labels = c(expression(paste(gamma, "00")), expression(paste(gamma, "10")), expression(paste(gamma, "20")), expression(paste(gamma, "30")), expression(paste(gamma, "40")), expression(paste(gamma, "50")), expression(paste(gamma, "60")), expression(paste(gamma, "70")), expression(paste(gamma, "01")), expression(paste(gamma, "02")), expression(paste(gamma, "11")), expression(paste(gamma, "21")), expression(paste(gamma, "32")), expression(paste(upsilon, "1")), expression(paste(upsilon, "2")), expression(paste(upsilon, "3")))) + 
    geom_vline(xintercept = 0, color = "gray40") +
    # geom_vline(xintercept = .10, linetype = "dashed", color = "gray40") +
    # geom_vline(xintercept = -.1, linetype = "dashed", color = "gray40") +
    geom_segment(y = 0, yend = 1.5, x = -.1, linetype = "dashed", color = "gray40") +
    geom_segment(y = 0, yend = 1.5, x = .1, linetype = "dashed", color = "gray40") +
    geom_segment(y = 1.5, yend = 8.5, x = -.05, linetype = "dashed", color = "gray40") +
    geom_segment(y = 1.5, yend = 8.5, x = .05, linetype = "dashed", color = "gray40") +
    geom_segment(y = 8.5, yend = 10.5, x = -.05, linetype = "dashed", color = "gray40") +
    geom_segment(y = 8.5, yend = 10.5, x = .05, linetype = "dashed", color = "gray40") +
    geom_segment(y = 10.5, yend = 13.5, x = -.035, linetype = "dashed", color = "gray40") +
    geom_segment(y = 10.5, yend = 13.5, x = .035, linetype = "dashed", color = "gray40") +
    geom_segment(y = 13.5, yend = 16.5, x = -.1, linetype = "dashed", color = "gray40") +
    geom_segment(y = 13.5, yend = 16.5, x = .1, linetype = "dashed", color = "gray40") +
    # scale_color_viridis_d() +
    # scale_color_brewer(palette = "Dark2") +
    # scale_color_manual(values = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#E6AB02", "#666666")) +
    # scale_color_lancet() +
    scale_color_manual(values = cbbPalette) +
    theme_minimal() +
    theme(panel.border = element_rect(colour = "gray25", fill = NA, size = .5), axis.text.x = element_text(size = 7, angle = 45, hjust = 1), axis.text.y = element_text(size = 9), legend.position = "bottom", legend.text = element_text(size = 9), axis.title = element_text(size = 12), legend.title = element_text(size = 12), panel.grid = element_line(color = "gray80")) +
    labs(
        x = "Bias",
        y = "Term",
        color = "Method"
    )
ggsave(bias1, file = "thesis/graphs/bias1.png",  width = 23, height = 30, units = "cm", dpi = "retina", bg = "white")
# eij and u0
bias2 <- bias_combined %>%
    filter(Term == "eij" | Term == "u0") %>%
    ggplot(aes(
        x = Bias,
        y = Term,
        color = Method
    )) +
    geom_point(position = position_jitter(seed = 123, width = 0, height = .2)) +
    geom_errorbar(aes(xmin = Bias - MCSE, xmax = Bias + MCSE), position = position_jitter(seed = 123, width = 0, height = .2), width = .25) +
    facet_grid(cols = vars(`Number of groups`, `Group size`), rows = vars(`Missingness mechanism`), labeller = labeller(.rows = label_value, .cols = label_both)) +
    scale_x_continuous(n.breaks = 15, minor_breaks = seq(-40, 60, 2.5)) +
    scale_y_discrete(limits = c("eij", "u0"), labels = c(expression(paste(epsilon, "ij")), expression(paste(upsilon, "0")))) +
    geom_vline(xintercept = 0, color = "gray40") +
    geom_segment(y = 0, yend = 1.5, x = -2.5, linetype = "dashed", color = "gray40") +
    geom_segment(y = 0, yend = 1.5, x = 2.5, linetype = "dashed", color = "gray40") +
    geom_segment(y = 1.5, yend = 3, x = -8.474903, linetype = "dashed", color = "gray40") +
    geom_segment(y = 1.5, yend = 3, x = 8.474903, linetype = "dashed", color = "gray40") +
    scale_color_manual(values = cbbPalette) +
    theme_minimal() +
    theme(panel.border = element_rect(colour = "gray25", fill = NA, size = .5), axis.text.x = element_text(size = 9, angle = 45, hjust = 1), axis.text.y = element_text(size = 9), legend.position = "bottom", legend.text = element_text(size = 9), axis.title = element_text(size = 12), legend.title = element_text(size = 12), panel.grid = element_line(color = "gray80")) +
    labs(
        x = "Absolute bias",
        y = "Term",
        color = "Method"
    )
ggsave(bias2, file = "thesis/graphs/bias2.png", width = 23, height = 9, units = "cm", dpi = "retina", bg = "white")
# Intercept 
biasintercept <- bias_combined %>%
    filter(Term == "beta0j") %>%
    ggplot(aes(
        x = Bias,
        y = Term,
        color = Method
    )) +
    geom_point(position = position_jitter(seed = 123, width = 0, height = .3)) +
    geom_errorbar(aes(xmin = Bias - MCSE, xmax = Bias + MCSE), position = position_jitter(seed = 123, width = 0, height = .3), width = .08) +
    facet_grid(cols = vars(`Number of groups`, `Group size`), rows = vars(`Missingness mechanism`), labeller = labeller(.rows = label_value, .cols = label_both)) +
    scale_x_continuous(n.breaks = 10, minor_breaks = seq(-.9, .7, .05)) +
    scale_y_discrete(limits = c("beta0j"), labels = c(expression(paste(gamma, "00")))) + 
    geom_vline(xintercept = 0, color = "gray40") +
    geom_segment(y = 0, yend = 2, x = -.1, linetype = "dashed", color = "gray40") +
    geom_segment(y = 0, yend = 2, x = .1, linetype = "dashed", color = "gray40") +
    scale_color_manual(values = cbbPalette) +
    theme_minimal() +
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
    geom_point(position = position_jitter(seed = 123, width = 0, height = .2)) +
    geom_errorbar(aes(xmin = Bias - MCSE, xmax = Bias + MCSE), position = position_jitter(seed = 123, width = 0, height = .2), width = .3) +
    facet_grid(cols = vars(`Number of groups`, `Group size`), rows = vars(`Missingness mechanism`), labeller = labeller(.rows = label_value, .cols = label_both)) +
    scale_x_continuous(n.breaks = 10) +
    scale_y_discrete(limits = c("beta1j", "beta2j", "beta3j", "beta4j", "beta5j", "beta6j", "beta7j"), labels = c(expression(paste(gamma, "10")), expression(paste(gamma, "20")), expression(paste(gamma, "30")), expression(paste(gamma, "40")), expression(paste(gamma, "50")), expression(paste(gamma, "60")), expression(paste(gamma, "70")))) + 
    geom_vline(xintercept = 0, color = "gray40") +
    geom_vline(xintercept = -.05, linetype = "dashed", color = "gray40") +
    geom_vline(xintercept = .05, linetype = "dashed", color = "gray40") +
    scale_color_manual(values = cbbPalette) +
    theme_minimal() +
    theme(panel.border = element_rect(colour = "gray25", fill = NA, size = .5), axis.text.x = element_text(size = 9, angle = 45, hjust = 1), axis.text.y = element_text(size = 9), legend.position = "bottom", legend.text = element_text(size = 9), axis.title = element_text(size = 12), legend.title = element_text(size = 12), panel.grid = element_line(color = "gray80")) +
    labs(
        x = "Absolute bias",
        y = "Term",
        color = "Method"
    )
ggsave(biaslevel1, file = "thesis/graphs/biaslevel1.png",  width = 23, height = 15, units = "cm", dpi = "retina", bg = "white")
# Fixed level-2 effects
biaslevel2 <- bias_combined %>%
    filter(Term == "z1" | Term == "z2") %>%
    ggplot(aes(
        x = Bias,
        y = Term,
        color = Method
    )) +
    geom_point(position = position_jitter(seed = 123, width = 0, height = .2)) +
    geom_errorbar(aes(xmin = Bias - MCSE, xmax = Bias + MCSE), position = position_jitter(seed = 123, width = 0, height = .2), width = .25) +
    facet_grid(cols = vars(`Number of groups`, `Group size`), rows = vars(`Missingness mechanism`), labeller = labeller(.rows = label_value, .cols = label_both)) +
    scale_x_continuous(n.breaks = 10, minor_breaks = seq(-.7, .9, .05)) +
    scale_y_discrete(limits = c("z1", "z2"), labels = c(expression(paste(gamma, "01")), expression(paste(gamma, "02")))) + 
    geom_vline(xintercept = 0, color = "gray40") +
    geom_segment(y = 0, yend = 3, x = -.05, linetype = "dashed", color = "gray40") +
    geom_segment(y = 0, yend = 3, x = .05, linetype = "dashed", color = "gray40") +
    scale_color_manual(values = cbbPalette) +
    theme_minimal() +
    theme(panel.border = element_rect(colour = "gray25", fill = NA, size = .5), axis.text.x = element_text(size = 9, angle = 45, hjust = 1), axis.text.y = element_text(size = 9), legend.position = "bottom", legend.text = element_text(size = 9), axis.title = element_text(size = 12), legend.title = element_text(size = 12), panel.grid = element_line(color = "gray80")) +
    labs(
        x = "Absolute bias",
        y = "Term",
        color = "Method"
    )
ggsave(biaslevel2, file = "thesis/graphs/biaslevel2.png",  width = 23, height = 9, units = "cm", dpi = "retina", bg = "white")
# Fixed cross-level effects
biascrosslevel <- bias_combined %>%
    filter(Term == "x1:z1" | Term == "x2:z1" | Term == "x3:z2") %>%
    ggplot(aes(
        x = Bias,
        y = Term,
        color = Method
    )) +
    geom_point(position = position_jitter(seed = 123, width = 0, height = .1)) +
    geom_errorbar(aes(xmin = Bias - MCSE, xmax = Bias + MCSE), position = position_jitter(seed = 123, width = 0, height = .1), width = .25) +
    facet_grid(cols = vars(`Number of groups`, `Group size`), rows = vars(`Missingness mechanism`), labeller = labeller(.rows = label_value, .cols = label_both)) +
    scale_x_continuous(n.breaks = 10, minor_breaks = seq(-.325, .1, .025)) +
    scale_y_discrete(limits = c("x1:z1", "x2:z1", "x3:z2"), labels = c(expression(paste(gamma, "11")), expression(paste(gamma, "21")), expression(paste(gamma, "32")))) + 
    geom_vline(xintercept = 0, color = "gray40") +
    geom_segment(y = 0, yend = 4, x = -.035, linetype = "dashed", color = "gray40") +
    geom_segment(y = 0, yend = 4, x = .035, linetype = "dashed", color = "gray40") +
    scale_color_manual(values = cbbPalette) +
    theme_minimal() +
    theme(panel.border = element_rect(colour = "gray25", fill = NA, size = .5), axis.text.x = element_text(size = 9, angle = 45, hjust = 1), axis.text.y = element_text(size = 9), legend.position = "bottom", legend.text = element_text(size = 9), axis.title = element_text(size = 12), legend.title = element_text(size = 12), panel.grid = element_line(color = "gray80")) +
    labs(
        x = "Absolute bias",
        y = "Term",
        color = "Method"
    )
ggsave(biascrosslevel, file = "thesis/graphs/biascrosslevel.png",  width = 23, height = 10, units = "cm", dpi = "retina", bg = "white")
# Random slopes
biasrandom <- bias_combined %>%
    filter(Term == "u1" | Term == "u2" | Term == "u3") %>%
    ggplot(aes(
        x = Bias,
        y = Term,
        color = Method
    )) +
    geom_point(position = position_jitter(seed = 123, width = 0, height = .1)) +
    geom_errorbar(aes(xmin = Bias - MCSE, xmax = Bias + MCSE), position = position_jitter(seed = 123, width = 0, height = .1), width = .25) +
    facet_grid(cols = vars(`Number of groups`, `Group size`), rows = vars(`Missingness mechanism`), labeller = labeller(.rows = label_value, .cols = label_both)) +
    scale_x_continuous(n.breaks = 10, minor_breaks = seq(-.6, .25, .025)) +
    scale_y_discrete(limits = c("u1", "u2", "u3"), labels = c(expression(paste(upsilon, "1")), expression(paste(upsilon, "2")), expression(paste(upsilon, "3")))) + 
    geom_vline(xintercept = 0, color = "gray40") +
    geom_segment(y = 0, yend = 4, x = -.1, linetype = "dashed", color = "gray40") +
    geom_segment(y = 0, yend = 4, x = .1, linetype = "dashed", color = "gray40") +
    scale_color_manual(values = cbbPalette) +
    theme_minimal() +
    theme(panel.border = element_rect(colour = "gray25", fill = NA, size = .5), axis.text.x = element_text(size = 9, angle = 45, hjust = 1), axis.text.y = element_text(size = 9), legend.position = "bottom", legend.text = element_text(size = 9), axis.title = element_text(size = 12), legend.title = element_text(size = 12), panel.grid = element_line(color = "gray80")) +
    labs(
        x = "Absolute bias",
        y = "Term",
        color = "Method"
    )
ggsave(biasrandom, file = "thesis/graphs/biasrandom.png",  width = 23, height = 10, units = "cm", dpi = "retina", bg = "white")
############
# Coverage #
############
coverage <- coverage_combined %>%
    ggplot(aes(
        x = Coverage,
        y = Term,
        color = Method
    )) +
    geom_point(position = position_jitter(seed = 123, width = 0, height = .1)) +
    geom_errorbar(aes(xmin = Coverage - MCSE, xmax = Coverage + MCSE), position = position_jitter(seed = 123, width = 0, height = .1), width = .3) +
    geom_vline(xintercept = .95, color = "gray40") +
    geom_vline(xintercept = .925, linetype = "dashed", color = "gray40") +
    geom_vline(xintercept = .975, linetype = "dashed", color = "gray40") +
    facet_grid(cols = vars(`Number of groups`, `Group size`), rows = vars(`Missingness mechanism`), labeller = labeller(.rows = label_value, .cols = label_both)) +
    scale_x_continuous(n.breaks = 10, minor_breaks = seq(0, 1, .025)) +
    scale_y_discrete(limits = c("beta0j", "beta1j", "beta2j", "beta3j", "beta4j", "beta5j", "beta6j", "beta7j", "z1", "z2", "x1:z1", "x2:z1", "x3:z2"), labels = c(expression(paste(gamma, "00")), expression(paste(gamma, "10")), expression(paste(gamma, "20")), expression(paste(gamma, "30")), expression(paste(gamma, "40")), expression(paste(gamma, "50")), expression(paste(gamma, "60")), expression(paste(gamma, "70")), expression(paste(gamma, "01")), expression(paste(gamma, "02")), expression(paste(gamma, "11")), expression(paste(gamma, "21")), expression(paste(gamma, "32")))) +
    scale_color_manual(values = cbbPalette) +
    theme_minimal() +
    theme(panel.border = element_rect(colour = "gray25", fill = NA, size = .5), axis.text.x = element_text(size = 9, angle = 45, hjust = 1), axis.text.y = element_text(size = 9), legend.position = "bottom", legend.text = element_text(size = 9), axis.title = element_text(size = 12), legend.title = element_text(size = 12), panel.grid = element_line(color = "gray80")) +
    labs(
        x = "Coverage",
        y = "Term",
        color = "Method"
    )
ggsave(coverage, file = "thesis/graphs/coverage.png", width = 23, height = 30, units = "cm", dpi = "retina", bg = "white")
# Fixed level-1 effects
coveragelevel1 <- coverage_combined %>%
    filter(Term == "beta0j" | Term == "beta1j" | Term == "beta2j" | Term == "beta3j" | Term == "beta4j" | Term == "beta5j" | Term == "beta6j" | Term == "beta7j") %>%
    ggplot(aes(
        x = Coverage,
        y = Term,
        color = Method
    )) +
    geom_point(position = position_jitter(seed = 123, width = 0, height = .2)) +
    geom_errorbar(aes(xmin = Coverage - MCSE, xmax = Coverage + MCSE), position = position_jitter(seed = 123, width = 0, height = .2), width = .3) +
    geom_vline(xintercept = .95, color = "gray40") +
    geom_vline(xintercept = .925, linetype = "dashed", color = "gray40") +
    geom_vline(xintercept = .975, linetype = "dashed", color = "gray40") +
    facet_grid(cols = vars(`Number of groups`, `Group size`), rows = vars(`Missingness mechanism`), labeller = labeller(.rows = label_value, .cols = label_both)) +
    scale_x_continuous(n.breaks = 10, minor_breaks = seq(0, 1, .025)) +
    scale_y_discrete(limits = c("beta0j", "beta1j", "beta2j", "beta3j", "beta4j", "beta5j", "beta6j", "beta7j"), labels = c(expression(paste(gamma, "00")), expression(paste(gamma, "10")), expression(paste(gamma, "20")), expression(paste(gamma, "30")), expression(paste(gamma, "40")), expression(paste(gamma, "50")), expression(paste(gamma, "60")), expression(paste(gamma, "70")))) +
    scale_color_manual(values = cbbPalette) +
    theme_minimal() +
    theme(panel.border = element_rect(colour = "gray25", fill = NA, size = .5), axis.text.x = element_text(size = 9, angle = 45, hjust = 1), axis.text.y = element_text(size = 9), legend.position = "bottom", legend.text = element_text(size = 9), axis.title = element_text(size = 12), legend.title = element_text(size = 12), panel.grid = element_line(color = "gray80")) +
    labs(
        x = "Coverage",
        y = "Term",
        color = "Method"
    )
ggsave(coveragelevel1, file = "thesis/graphs/coveragelevel1.png", width = 23, height = 15, units = "cm", dpi = "retina", bg = "white")
# Fixed level-2 effects
coveragelevel2 <- coverage_combined %>%
    filter(Term == "z1" | Term == "z2" | Term == "x1:z1" | Term == "x2:z1" | Term == "x3:z2") %>%
    ggplot(aes(
        x = Coverage,
        y = Term,
        color = Method
    )) +
    geom_point(position = position_jitter(seed = 123, width = 0, height = .2)) +
    geom_errorbar(aes(xmin = Coverage - MCSE, xmax = Coverage + MCSE), position = position_jitter(seed = 123, width = 0, height = .2), width = .25) +
    geom_vline(xintercept = .95, color = "gray40") +
    geom_vline(xintercept = .925, linetype = "dashed", color = "gray40") +
    geom_vline(xintercept = .975, linetype = "dashed", color = "gray40") +
    facet_grid(cols = vars(`Number of groups`, `Group size`), rows = vars(`Missingness mechanism`), labeller = labeller(.rows = label_value, .cols = label_both)) +
    scale_x_continuous(n.breaks = 10, minor_breaks = seq(0, 1, .025)) +
    scale_y_discrete(limits = c("z1", "z2", "x1:z1", "x2:z1", "x3:z2"), labels = c(expression(paste(gamma, "01")), expression(paste(gamma, "02")), expression(paste(gamma, "11")), expression(paste(gamma, "21")), expression(paste(gamma, "32")))) +
    scale_color_manual(values = cbbPalette) +
    theme_minimal() +
    theme(panel.border = element_rect(colour = "gray25", fill = NA, size = .5), axis.text.x = element_text(size = 9, angle = 45, hjust = 1), axis.text.y = element_text(size = 9), legend.position = "bottom", legend.text = element_text(size = 9), axis.title = element_text(size = 12), legend.title = element_text(size = 12), panel.grid = element_line(color = "gray80")) +
    labs(
        x = "Coverage",
        y = "Term",
        color = "Method"
    )
ggsave(coveragelevel2, file = "thesis/graphs/coveragelevel2.png", width = 23, height = 13, units = "cm", dpi = "retina", bg = "white")
#######
# CIW #
#######
ciw <- ciw_combined %>%
    ggplot(aes(
        x = CIW,
        y = Term,
        color = Method
        # fill = Method
    )) +
    geom_point(position = position_jitter(seed = 123, width = 0, height = .1)) +
    facet_grid(cols = vars(`Number of groups`, `Group size`), rows = vars(`Missingness mechanism`), labeller = labeller(.rows = label_value, .cols = label_both)) +
    scale_x_continuous(n.breaks = 10, minor_breaks = seq(0, 8, .25)) +
    scale_y_discrete(limits = c("beta0j", "beta1j", "beta2j", "beta3j", "beta4j", "beta5j", "beta6j", "beta7j", "z1", "z2", "x1:z1", "x2:z1", "x3:z2"), labels = c(expression(paste(gamma, "00")), expression(paste(gamma, "10")), expression(paste(gamma, "20")), expression(paste(gamma, "30")), expression(paste(gamma, "40")), expression(paste(gamma, "50")), expression(paste(gamma, "60")), expression(paste(gamma, "70")), expression(paste(gamma, "01")), expression(paste(gamma, "02")), expression(paste(gamma, "11")), expression(paste(gamma, "21")), expression(paste(gamma, "32")))) +
    scale_color_manual(values = cbbPalette) +
    theme_minimal() +
    theme(panel.border = element_rect(colour = "gray25", fill = NA, size = .5), axis.text.x = element_text(size = 9, angle = 0, hjust = 1), axis.text.y = element_text(size = 9), legend.position = "bottom", legend.text = element_text(size = 9), axis.title = element_text(size = 12), legend.title = element_text(size = 12), panel.grid = element_line(color = "gray80")) +
    labs(
        x = "CIW",
        y = "Term",
        color = "Method"
        # fill = "Method"
    )
ggsave(ciw, file = "thesis/graphs/ciw.png", width = 23, height = 30, units = "cm", dpi = "retina", bg = "white")
# Intercept
ciwintercept <- ciw_combined %>%
    filter(Term == "beta0j") %>%
    ggplot(aes(
        x = CIW,
        y = Term,
        color = Method
        # fill = Method
    )) +
    geom_point(position = position_jitter(seed = 123, width = 0, height = .1)) +
    facet_grid(cols = vars(`Number of groups`, `Group size`), rows = vars(`Missingness mechanism`), labeller = labeller(.rows = label_value, .cols = label_both)) +
    scale_x_continuous(n.breaks = 5, minor_breaks = seq(4, 8, .25)) +
    scale_y_discrete(limits = c("beta0j"), labels = c(expression(paste(gamma, "00")))) +
    scale_color_manual(values = cbbPalette) +
    theme_minimal() +
    theme(panel.border = element_rect(colour = "gray25", fill = NA, size = .5), axis.text.x = element_text(size = 9, angle = 0, hjust = 1), axis.text.y = element_text(size = 9), legend.position = "bottom", legend.text = element_text(size = 9), axis.title = element_text(size = 12), legend.title = element_text(size = 12), panel.grid = element_line(color = "gray80")) +
    labs(
        x = "CIW",
        y = "Term",
        color = "Method"
    )
ggsave(ciwintercept, file = "thesis/graphs/ciwintercept.png", width = 23, height = 9, units = "cm", dpi = "retina", bg = "white")
# Fixed level-1 effects
ciwlevel1 <- ciw_combined %>%
    filter(Term == "beta1j" | Term == "beta2j" | Term == "beta3j" | Term == "beta4j" | Term == "beta5j" | Term == "beta6j" | Term == "beta7j") %>%
    ggplot(aes(
        x = CIW,
        y = Term,
        color = Method
        # fill = Method
    )) +
    geom_point(position = position_jitter(seed = 123, width = 0, height = .1)) +
    facet_grid(cols = vars(`Number of groups`, `Group size`), rows = vars(`Missingness mechanism`), labeller = "label_both") +
    scale_x_continuous(n.breaks = 10) +
    scale_y_discrete(limits = c("beta1j", "beta2j", "beta3j", "beta4j", "beta5j", "beta6j", "beta7j"), labels = c(expression(paste(gamma, "10")), expression(paste(gamma, "20")), expression(paste(gamma, "30")), expression(paste(gamma, "40")), expression(paste(gamma, "50")), expression(paste(gamma, "60")), expression(paste(gamma, "70")))) +
    scale_color_manual(values = cbbPalette) +
    theme_minimal() +
    theme(panel.border = element_rect(colour = "gray25", fill = NA, size = .5), axis.text.x = element_text(size = 9, angle = 45, hjust = 1), axis.text.y = element_text(size = 9), legend.position = "bottom", legend.text = element_text(size = 9), axis.title = element_text(size = 12), legend.title = element_text(size = 12), panel.grid = element_line(color = "gray80")) +
    labs(
        x = "CIW",
        y = "Term",
        color = "Method"
    )
ggsave(ciwlevel1, file = "thesis/graphs/ciwlevel1.png", width = 23, height = 15, units = "cm", dpi = "retina", bg = "white")
# Fixed level-2 effects
ciwlevel2 <- ciw_combined %>%
    filter(Term == "z1" | Term == "z2") %>%
    ggplot(aes(
        x = CIW,
        y = Term,
        color = Method
        # fill = Method
    )) +
    geom_point(position = position_jitter(seed = 123, width = 0, height = .2)) +
    facet_grid(cols = vars(`Number of groups`, `Group size`), rows = vars(`Missingness mechanism`), labeller = labeller(.rows = label_value, .cols = label_both)) +
    scale_x_continuous(n.breaks = 10, minor_breaks = seq(0, 8, .25)) +
    scale_y_discrete(limits = c("z1", "z2"), labels = c(expression(paste(gamma, "01")), expression(paste(gamma, "02")))) +
    scale_color_manual(values = cbbPalette) +
    theme_minimal() +
    theme(panel.border = element_rect(colour = "gray25", fill = NA, size = .5), axis.text.x = element_text(size = 9, angle = 0, hjust = 1), axis.text.y = element_text(size = 9), legend.position = "bottom", legend.text = element_text(size = 9), axis.title = element_text(size = 12), legend.title = element_text(size = 12), panel.grid = element_line(color = "gray80")) +
    labs(
        x = "CIW",
        y = "Term",
        color = "Method"
    )
ggsave(ciwlevel2, file = "thesis/graphs/ciwlevel2.png", width = 23, height = 9, units = "cm", dpi = "retina", bg = "white")
# Cross-level effects
ciwcrosslevel <- ciw_combined %>%
    filter(Term == "x1:z1" | Term == "x2:z1" | Term == "x3:z2") %>%
    ggplot(aes(
        x = CIW,
        y = Term,
        color = Method
        # fill = Method
    )) +
    geom_point(position = position_jitter(seed = 123, width = 0, height = .2)) +
    facet_grid(cols = vars(`Number of groups`, `Group size`), rows = vars(`Missingness mechanism`), labeller = labeller(.rows = label_value, .cols = label_both)) +
    scale_x_continuous(n.breaks = 10) +
    scale_y_discrete(limits = c("x1:z1", "x2:z1", "x3:z2"), labels = c(expression(paste(gamma, "11")), expression(paste(gamma, "21")), expression(paste(gamma, "32")))) +
    scale_color_manual(values = cbbPalette) +
    theme_minimal() +
    theme(panel.border = element_rect(colour = "gray25", fill = NA, size = .5), axis.text.x = element_text(size = 9, angle = 45, hjust = 1), axis.text.y = element_text(size = 9), legend.position = "bottom", legend.text = element_text(size = 9), axis.title = element_text(size = 12), legend.title = element_text(size = 12), panel.grid = element_line(color = "gray80")) +
    labs(
        x = "CIW",
        y = "Term",
        color = "Method"
    )
ggsave(ciwcrosslevel, file = "thesis/graphs/ciwcrosslevel.png", width = 23, height = 10, units = "cm", dpi = "retina", bg = "white")
