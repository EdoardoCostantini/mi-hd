# Title:    Initialization scirpt, general (functions and packages)
# Project:  Imputing High Dimensional Data
# Author:   Edoardo Costantini
# Created:  2020-07-01
# Modified: 2023-03-20

# Packages ----------------------------------------------------------------

pack_list <- c("tidyverse",
            #    "CVTuningCov",
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
            #  "devtools", #TODO: DO I need this? creates problems with IVEware
               "rstan",
               "gridExtra",
               "grid",
               "plyr", # for round_any
               "plot.matrix", # for missing data pattern plot
               "xtable",
               "PcAux",
               "miceadds",
               "ridge",
               "imputeR",
               "blasso")

# If not available, install
for(PKG in pack_list){
    if (!require(PKG, character.only = TRUE)) install.packages(PKG)
}

# then load all
lapply(pack_list, library, character.only = TRUE)

# Support Functions -------------------------------------------------------

source("./functions.R")
source("./fun_DURR_impute.R")
source("./fun_IURR_impute.R")
source("./fun_blasso_impute.R")
source("./fun_bridge_impute.R")
source("./fun_PCA_impute.R")
source("./fun_IVEware.R")
source("./fun_MICE_qp.R")
source("./fun_MICE_TR_impute.R")
source("./fun_MICE_cp.R")
source("./fun_RANF_impute.R")
source("./fun_CART_impute.R")
source("./fun_missFor_impute.R")
source("./functions_genDt.R")
source("./subroutines.R")
source("./simMissingness.R")
