### Title:    Conergence Check Blasso
### Project:  Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-08-24

rm(list = ls())
source("./init_general.R")

# Results -----------------------------------------------------------------
# Load data
dataDir <- "../output/"

filename <- "exp4_cv_bridge_20201015_1558" # k btw 1e4 and 1e-6
filename <- "exp4_cv_bridge_20201015_2037" # k btw 1e-1 and 1e-8
filename <- "exp4_cv_bridge_20201027_1411" # k btw 500 and 1500

out <- readRDS(paste0(dataDir, filename, ".rds")) 

# Checked Range
range(out$conds$ridge)

# Obtain conditions with cv ridge
conds_bridge <- bridge_cv(out)
conds_bridge

# Comparisons -------------------------------------------------------------
# Narrowing down areound condition 2 ridge penalty
# by using values around 1e3 to check best ridge penalty there.
filename <- "exp4_cv_bridge_20220223_1646" # k btw 1e-1 and 1e-8

out <- readRDS(paste0(dataDir, filename, ".rds")) 

# Checked Range
range(out$conds$ridge)

# Obtain conditions with cv ridge
conds_bridge <- bridge_cv(out, exp_factors = c("n"))
conds_bridge$values
conds_bridge$plot
