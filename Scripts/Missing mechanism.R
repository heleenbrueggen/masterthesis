#####################
# Missing mechanism #
#####################
#############
# Libraries #
#############
library(tidyverse)
library(mice)
library(readr)
library(magrittr)
################
# Setting seed #
################
set.seed(123)
################
# Setting path #
################
path <- "/Volumes/Heleen 480GB/MBART-MICE files/"
###########################
# Defining design factors #
###########################
combinations <- read_rds(paste(path, "data/combinations.rds", sep = ""))
#############################
# Storing names of datasets #
#############################
names <- read_rds(paste(path, "data/names.rds", sep = ""))
#############
# Load data #
#############
# Extracting amputation objects from results
miss <- list()
for (i in seq_len(nrow(combinations))) {
    miss[[i]] <- read_rds(paste(path, "data/missing/miss_", names[i], ".rds", sep = "")) %>% map(~ .x %$% miss)
}
###########
# bwplots #
###########
bwplot(miss[[8]][[1]])
