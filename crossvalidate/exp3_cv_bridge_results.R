### Title:    Conergence Check Blasso
### Project:  Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-08-24

rm(list = ls())
source("./init_general.R")

# Load data
out <- readRDS("../output/exp3_cv_bridge_20200824_1618.rds")

# Obtain conditions with cv ridge
conds_bridge <- bridge_cv(out)
conds_bridge