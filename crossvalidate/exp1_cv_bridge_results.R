# Title:    Simulation script for ridge paramter cross-validation
# Project:  Imputing High Dimensional Data
# Author:   Edoardo Costantini
# Created:  2020-08-25
# Modified: 2023-04-08

rm(list = ls())
source("./init_general.R")

# Load data
out <- readRDS("../output/exp1_cv_bridge_20220224_1042.rds")

# Obtain conditions with cv ridge
conds_bridge <- cvParm(
    out,
    cv.parm = "ridge",
    exp_factors = colnames(out$conds)[c(2, 4)]
)
conds_bridge$solution
conds_bridge$solution_1se
conds_bridge$plot
