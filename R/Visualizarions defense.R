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
##############
# Formatting #
##############
bias_combined <- bias_combined %>%
    select(Term, Bias, MCSE, Method) %>%
    group_by(Term, Method) %>%
    reframe(
        MCSE = mean(MCSE),
        Bias = mean(Bias)
    )
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
    geom_point(position = position_jitter(seed = 123, width = 0, height = .2), size = 2) +
    geom_errorbar(aes(xmin = Bias - MCSE, xmax = Bias + MCSE), position = position_jitter(seed = 123, width = 0, height = .2), width = .1) +
    # facet_grid(cols = vars(`Number of groups`, `Group size`), rows = vars(`Missingness mechanism`), labeller = labeller(.rows = label_value, .cols = label_both)) +
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
ggsave(bias2, file = "defensepresentation/graphs/bias2.png", width = 23, height = 9, units = "cm", dpi = "retina", bg = NULL)
# Fixed level-2 effects
biaslevel2 <- bias_combined %>%
    filter(Term == "z1" | Term == "z2") %>%
    ggplot(aes(
        x = Bias,
        y = Term,
        color = Method
    )) +
    geom_point(position = position_jitter(seed = 123, width = 0, height = .2), size = 2) +
    geom_errorbar(aes(xmin = Bias - MCSE, xmax = Bias + MCSE), position = position_jitter(seed = 123, width = 0, height = .2), width = .1) +
    # facet_grid(cols = vars(`Number of groups`, `Group size`), rows = vars(`Missingness mechanism`), labeller = labeller(.rows = label_value, .cols = label_both)) +
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
ggsave(biaslevel2, file = "defensepresentation/graphs/biaslevel2.png",  width = 23, height = 9, units = "cm", dpi = "retina", bg = NULL)
# Fixed cross-level effects
biascrosslevel <- bias_combined %>%
    filter(Term == "x1:z1" | Term == "x2:z1" | Term == "x3:z2") %>%
    ggplot(aes(
        x = Bias,
        y = Term,
        color = Method
    )) +
    geom_point(position = position_jitter(seed = 123, width = 0, height = .1), size = 2) +
    geom_errorbar(aes(xmin = Bias - MCSE, xmax = Bias + MCSE), position = position_jitter(seed = 123, width = 0, height = .1), width = .1) +
    # facet_grid(cols = vars(`Number of groups`, `Group size`), rows = vars(`Missingness mechanism`), labeller = labeller(.rows = label_value, .cols = label_both)) +
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
ggsave(biascrosslevel, file = "defensepresentation/graphs/biascrosslevel.png",  width = 23, height = 10, units = "cm", dpi = "retina", bg = NULL)
# Random slopes
biasrandom <- bias_combined %>%
    filter(Term == "u1" | Term == "u2" | Term == "u3") %>%
    ggplot(aes(
        x = Bias,
        y = Term,
        color = Method
    )) +
    geom_point(position = position_jitter(seed = 123, width = 0, height = .1), size = 2) +
    geom_errorbar(aes(xmin = Bias - MCSE, xmax = Bias + MCSE), position = position_jitter(seed = 123, width = 0, height = .1), width = .1) +
    # facet_grid(cols = vars(`Number of groups`, `Group size`), rows = vars(`Missingness mechanism`), labeller = labeller(.rows = label_value, .cols = label_both)) +
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
ggsave(biasrandom, file = "defensepresentation/graphs/biasrandom.png",  width = 23, height = 10, units = "cm", dpi = "retina", bg = NULL)
