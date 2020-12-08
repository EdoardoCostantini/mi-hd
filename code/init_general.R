### Title:    Initialization scirpt, general (functions and packages)
### Project:  Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-07-01

# Packages ----------------------------------------------------------------

pack_list <- c("tidyverse",
               "CVTuningCov",
               "mvtnorm",
               "monomvn",
               "glmnet",
               "rlecuyer",
               "parallel",
               "caret", 
               "missForest",
               "truncnorm",
               "lavaan",
               "FactoMineR",
               "devtools",
               "xtable",
               "rstan",
               "gridExtra",
               "grid",
               "plyr", # for round_any
               "xtable",
               "PcAux",
               "blasso")

lapply(pack_list, library, character.only = TRUE)

# library(tidyverse)
# library(CVTuningCov) # for AR1 function
# library(mvtnorm)
# library(monomvn)
# library(glmnet)
# library(rlecuyer) # for seed in parallel
# library(parallel)
# library(caret) # for Elastic cv functions
# library(missForest) # for missForest imputation
# library(truncnorm) # for truncnorm sampling in probit models
# library(lavaan)
# library(FactoMineR)
# library(devtools)
# library(xtable)
# library(rstan) # for Rhat computation
# 
# # Special pacakges
# library(PcAux)  # requires installation with devtools
# library(blasso) # requires installation with zipped folder

# Support Functions -------------------------------------------------------

source("./functions.R")
source("./functions_PCA_impute.R")
source("./fun_DURR_impute.R")
source("./fun_IURR_impute.R")
source("./fun_blasso_impute.R")
source("./fun_bridge_impute.R")
source("./fun_PCA_impute.R")
source("./fun_MICE_CART_impute.R")
source("./fun_MICE_RF_impute.R")
source("./fun_MICE_TR_impute.R")
source("./fun_RANF_impute.R")
source("./fun_CART_impute.R")
source("./fun_missFor_impute.R")
source("./functions_genDt.R")
source("./subroutines.R")
source("./simMissingness.R")
source("./functions_EVS.R")

# Experiment Specific -----------------------------------------------------
# source("./gen_lavaan_model.R") # generate txt file for lavaan model
