#################################
# Formatting for visualizations #
#################################
#############
# Libraries #
#############
library(tidyverse)
################
# Setting seed #
################
set.seed(123)
################
# Setting path # 
################
path <- "/Volumes/Heleen 480GB/MBART-MICE files/"
############################
# Functions for formatting # 
############################
format.bias <- function(x, y, method, name = c("Bias", "MCSE")) {
    data <- cbind(x, y) %>% 
        pivot_longer(cols = beta0j:eij, names_to = "term", values_to = name) %>%
        mutate(ngroup = as.factor(ngroup), groupsize = as.factor(groupsize), icc = as.factor(icc), mar_mcar = as.factor(mar_mcar), miss = as.factor(miss), g = as.factor(g), method = method)
    colnames(data) <- c("Number of groups", "Group size", "ICC", "Missingness mechanism", "Percentage of missing data", "Gamma", "Term", name, "Method")
    data <- data %>%
        mutate(`Missingness mechanism` = ifelse(`Missingness mechanism` == "mar", "MAR", "MCAR"))
    return(data)
}
format.coverage <- function(x, y, method, name = c("Coverage", "MCSE")) {
    data <- cbind(x, y) %>%
        pivot_longer(cols = beta0j:`x3:z2`, names_to = "term", values_to = name) %>%
        mutate(ngroup = as.factor(ngroup), groupsize = as.factor(groupsize), icc = as.factor(icc), mar_mcar = as.factor(mar_mcar), miss = as.factor(miss), g = as.factor(g), method = method)
    colnames(data) <- c("Number of groups", "Group size", "ICC", "Missingness mechanism", "Percentage of missing data", "Gamma", "Term", name, "Method")
    data <- data %>%
        mutate(`Missingness mechanism` = ifelse(`Missingness mechanism` == "mar", "MAR", "MCAR"))
    return(data)
}
format.ciw <- function(x, y, method, name = "CIW") {
    data <- cbind(x, y) %>%
        pivot_longer(cols = beta0j:`x3:z2`, names_to = "term", values_to = name) %>%
        mutate(ngroup = as.factor(ngroup), groupsize = as.factor(groupsize), icc = as.factor(icc), mar_mcar = as.factor(mar_mcar), miss = as.factor(miss), g = as.factor(g), method = method)
    colnames(data) <- c("Number of groups", "Group size", "ICC", "Missingness mechanism", "Percentage of missing data", "Gamma", "Term", name, "Method")
    data <- data %>%
        mutate(`Missingness mechanism` = ifelse(`Missingness mechanism` == "mar", "MAR", "MCAR"))
    return(data)
}
#####################
# Complete analysis #
#####################
#######################
# Defining parameters #
#######################
combinations <- read_rds(paste(path, "data/nomissing/combinations.rds", sep = ""))
# Bias
bias_nomiss <- read_rds(paste(path, "results/evaluations/bias_nomiss.rds", sep = "")) %>% 
    format.bias(., combinations, method = "true", name = "Bias") %>%
    cbind(read_rds(paste(path, "results/evaluations/mcse_bias_nomiss.rds", sep = "")) %>%
        format.bias(., combinations, method = "true", name = "MCSE") %>%
        dplyr::select(MCSE))
# Coverage
coverage_nomiss <- read_rds(paste(path, "results/evaluations/coverage_nomiss.rds", sep = "")) %>%
    format.coverage(., combinations, method = "true", name = "Coverage") %>%
    cbind(read_rds(paste(path, "results/evaluations/mcse_coverage_nomiss.rds", sep = "")) %>%
        format.coverage(., combinations, method = "true", name = "MCSE") %>%
        dplyr::select(MCSE))
# CIW 
ciw_nomiss <- read_rds(paste(path, "results/evaluations/ciw_nomiss.rds", sep = "")) %>%
    format.ciw(., combinations, method = "true", name = "CIW") 
#############################################
# Defining parameters for data with missing #
#############################################
combinations.imp <- read_rds(paste(path, "data/combinations.rds", sep = ""))
#####################
# Listwise deletion #
#####################
# Bias
bias_ld <- read_rds(paste(path, "results/evaluations/bias_ld.rds", sep = "")) %>%
    format.bias(., combinations.imp, method = "ld", name = "Bias") %>%
    cbind(read_rds(paste(path, "results/evaluations/mcse_bias_ld.rds", sep = "")) %>% 
    format.bias(., combinations.imp, method = "ld", name = "MCSE") %>% 
    dplyr::select(MCSE))
# Coverage
coverage_ld <- read_rds(paste(path, "results/evaluations/coverage_ld.rds", sep = "")) %>%
    format.coverage(., combinations.imp, method = "ld", name = "Coverage") %>%
    cbind(read_rds(paste(path, "results/evaluations/mcse_coverage_ld.rds", sep = "")) %>% 
    format.coverage(., combinations.imp, method = "ld", name = "MCSE") %>% 
    dplyr::select(MCSE))
# CIW 
ciw_ld <- read_rds(paste(path, "results/evaluations/ciw_ld.rds", sep = "")) %>%
    format.ciw(., combinations.imp, method = "ld", name = "CIW")
#######
# PMM # 
#######
# Bias 
bias_pmm <- read_rds(paste(path, "results/evaluations/bias_pmm.rds", sep = "")) %>%
    format.bias(., combinations.imp, method = "pmm", name = "Bias") %>%
    cbind(read_rds(paste(path, "results/evaluations/mcse_bias_pmm.rds", sep = "")) %>% 
    format.bias(., combinations.imp, method = "pmm", name = "MCSE") %>% 
    dplyr::select(MCSE))
# Coverage
coverage_pmm <- read_rds(paste(path, "results/evaluations/coverage_pmm.rds", sep = "")) %>%
    format.coverage(., combinations.imp, method = "pmm", name = "Coverage") %>%
    cbind(read_rds(paste(path, "results/evaluations/mcse_coverage_pmm.rds", sep = "")) %>% 
    format.coverage(., combinations.imp, method = "pmm", name = "MCSE") %>% 
    dplyr::select(MCSE))
# CIW 
ciw_pmm <- read_rds(paste(path, "results/evaluations/ciw_pmm.rds", sep = "")) %>%
    format.ciw(., combinations.imp, method = "pmm", name = "CIW")
##########
# 2l.PMM # 
##########
# Bias
bias_2l.pmm <- read_rds(paste(path, "results/evaluations/bias_2l.pmm.rds", sep = "")) %>%
    format.bias(., combinations.imp, method = "2l.pmm", name = "Bias") %>%
    cbind(read_rds(paste(path, "results/evaluations/mcse_bias_2l.pmm.rds", sep = "")) %>% 
    format.bias(., combinations.imp, method = "2l.pmm", name = "MCSE") %>% 
    dplyr::select(MCSE))
# Coverage
coverage_2l.pmm <- read_rds(paste(path, "results/evaluations/coverage_2l.pmm.rds", sep = "")) %>%
    format.coverage(., combinations.imp, method = "2l.pmm", name = "Coverage") %>%
    cbind(read_rds(paste(path, "results/evaluations/mcse_coverage_2l.pmm.rds", sep = "")) %>% 
    format.coverage(., combinations.imp, method = "2l.pmm", name = "MCSE") %>% 
    dplyr::select(MCSE))
# CIW
ciw_2l.pmm <- read_rds(paste(path, "results/evaluations/ciw_2l.pmm.rds", sep = "")) %>%
    format.ciw(., combinations.imp, method = "2l.pmm", name = "CIW")
########
# bart #
########
# Bias
bias_bart <- read_rds(paste(path, "results/evaluations/bias_bart.rds", sep = "")) %>%
    format.bias(., combinations.imp, method = "bart", name = "Bias") %>%
    cbind(read_rds(paste(path, "results/evaluations/mcse_bias_bart.rds", sep = "")) %>% 
    format.bias(., combinations.imp, method = "bart", name = "MCSE") %>% 
    dplyr::select(MCSE))
# Coverage
coverage_bart <- read_rds(paste(path, "results/evaluations/coverage_bart.rds", sep = "")) %>%
    format.coverage(., combinations.imp, method = "bart", name = "Coverage") %>%
    cbind(read_rds(paste(path, "results/evaluations/mcse_coverage_bart.rds", sep = "")) %>% 
    format.coverage(., combinations.imp, method = "bart", name = "MCSE") %>% 
    dplyr::select(MCSE))
# CIW
ciw_bart <- read_rds(paste(path, "results/evaluations/ciw_bart.rds", sep = "")) %>%
    format.ciw(., combinations.imp, method = "bart", name = "CIW")
#########
# rbart #
#########
# Bias
bias_rbart <- read_rds(paste(path, "results/evaluations/bias_rbart.rds", sep = "")) %>%
    format.bias(., combinations.imp, method = "rbart", name = "Bias") %>%
    cbind(read_rds(paste(path, "results/evaluations/mcse_bias_rbart.rds", sep = "")) %>% 
    format.bias(., combinations.imp, method = "rbart", name = "MCSE") %>% 
    dplyr::select(MCSE))
# Coverage
coverage_rbart <- read_rds(paste(path, "results/evaluations/coverage_rbart.rds", sep = "")) %>%
    format.coverage(., combinations.imp, method = "rbart", name = "Coverage") %>%
    cbind(read_rds(paste(path, "results/evaluations/mcse_coverage_rbart.rds", sep = "")) %>% 
    format.coverage(., combinations.imp, method = "rbart", name = "MCSE") %>% 
    dplyr::select(MCSE))
# CIW 
ciw_rbart <- read_rds(paste(path, "results/evaluations/ciw_rbart.rds", sep = "")) %>%
    format.ciw(., combinations.imp, method = "rbart", name = "CIW")
#############
# stan4bart #
#############
# Bias 
bias_stan4bart <- read_rds(paste(path, "results/evaluations/bias_stan4bart.rds", sep = "")) %>%
    format.bias(., combinations.imp, method = "stan4bart", name = "Bias") %>%
    cbind(read_rds(paste(path, "results/evaluations/mcse_bias_stan4bart.rds", sep = "")) %>% 
    format.bias(., combinations.imp, method = "stan4bart", name = "MCSE") %>% 
    dplyr::select(MCSE))
# Coverage
coverage_stan4bart <- read_rds(paste(path, "results/evaluations/coverage_stan4bart.rds", sep = "")) %>%
    format.coverage(., combinations.imp, method = "stan4bart", name = "Coverage") %>%
    cbind(read_rds(paste(path, "results/evaluations/mcse_coverage_stan4bart.rds", sep = "")) %>% 
    format.coverage(., combinations.imp, method = "stan4bart", name = "MCSE") %>% 
    dplyr::select(MCSE))
# CIW
ciw_stan4bart <- read_rds(paste(path, "results/evaluations/ciw_stan4bart.rds", sep = "")) %>%
    format.ciw(., combinations.imp, method = "stan4bart", name = "CIW")
#####################
# Combining results # 
#####################
bias_combined <- rbind(bias_ld, bias_nomiss, bias_pmm, bias_2l.pmm, bias_bart, bias_rbart, bias_stan4bart)
coverage_combined <- rbind(coverage_ld, coverage_nomiss, coverage_pmm, coverage_2l.pmm, coverage_bart, coverage_rbart, coverage_stan4bart)
ciw_combined <- rbind(ciw_ld, ciw_nomiss, ciw_pmm, ciw_2l.pmm, ciw_bart, ciw_rbart, ciw_stan4bart)
