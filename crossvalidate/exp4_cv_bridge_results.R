# Title:    Conergence Check Blasso
# Project:  Imputing High Dimensional Data
# Author:   Edoardo Costantini
# Created:  2020-08-24
# Modified:  2023-04-08

rm(list = ls())
source("./init_general.R")

# Results -----------------------------------------------------------------
# Load data
dataDir <- "../output/"
filename <- "exp4_cv_bridge_20220223_1646"
out <- readRDS(paste0(dataDir, filename, ".rds")) 

# Checked Range of rdige
range(out$conds$ridge)

# Obtain conditions with cv ridge
conds_bridge <- cvParm(out, exp_factors = c("n"))
conds_bridge$solution
conds_bridge$solution_1se
conds_bridge$plot
