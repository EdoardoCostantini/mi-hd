# Title:    Simulation script for ridge paramter cross-validation
# Project:  Imputing High Dimensional Data
# Author:   Anonymized for peer review
# Created:  2020-08-25
# Modified: 2022-02-25

rm(list = ls())
source("./init_general.R")

# Load data
out <- readRDS("../output/exp1_cv_bridge_20220224_1042.rds")

# Obtain conditions with cv ridge
conds_bridge <- bridge_cv(out,
                          exp_factors = colnames(out$conds)[c(2, 4)])
conds_bridge$values
conds_bridge$plot
