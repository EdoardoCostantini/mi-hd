### Title:    Conergence Check Blasso
### Project:  Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-08-24

rm(list = ls())
source("./init_general.R")

# Load data
dataDir <- "../output/"

filename <- "exp4_cv_bridge_20201015_1558" # k btw 1e4 and 1e-6
filename <- "exp4_cv_bridge_20201015_2037" # k btw 1e-1 and 1e-8

out <- readRDS(paste0(dataDir, filename, ".rds")) 

# Range
range(out$conds$ridge)

# Obtain conditions with cv ridge
conds_bridge <- bridge_cv(out)
conds_bridge
