### Title:    Replication Zhao Long 2016 - initialization script
### Author:   Edoardo Costantini
### Created:  2020-02-20
### Modified: 2020-02-20

# Packages ----------------------------------------------------------------

library(tidyverse)
library(CVTuningCov) 
library(mvtnorm)
library(monomvn)
library(glmnet)
library(rlecuyer) # for seed in parallel
library(MIBRR) # for bl function fitting blasso

# Support Functions -------------------------------------------------------

source("./functions.R")
source("./fun_DURR_impute.R")
source("./fun_IURR_impute.R")
source("./fun_BLasso_impute.R")
source("./fun_blassoHans_impute.R")
source("./subroutines.R")

# Fixed Parameters --------------------------------------------------------

parms <- list()
parms$dt_rep  <- 500 # 500 on original # replications for averaging results
parms$chains  <- 10  # number of imputed datasets to pool esitmaes from
parms$iters   <- 1   # number of iterations for the MI procedure
parms$n       <- 100 # number of cases for data generation
parms$y_var   <- 3   # error variance for y generation
parms$b       <- 1   # betas and intercept in linear model to gen y
parms$b_miss_model <- c(-1, -.1, 2, -2) # betas for generation of missingness
parms$stnr    <- c(1, sqrt(.2), sqrt(.08)) # signal to noise ratio
parms$S_all   <- list(q4=(c(2,3,50,51)-1),
                      q20=(c(2:11, 50:59)-1),
                      q50=(c(2:11, 50:59, 70:79, 90:99, 110:119)-1))
parms$formula <- "y ~ z1 + z2 + z3"
parms$alphaCI <- .95 # confidence level for parameters CI
parms$k       <- 3 # number of predictors for analysis model
parms$iter_bl <- 1e4 # 1e3 but actually no differences
parms$burn_bl <- 1 * parms$iter_bl

parms$k_IURR <- 0 # k value to bias coef sampling covariance matrix in IURR
                  # procedure to solve possible issues of singularity

parms$methods <- c("DURR", "IURR", "Blas_HS", "Blas_PC", "MI_T", "MI_50", "CC") 
  # for naming objects, hence order is important (CC always last)

# Replicability related
parms$seed     <- 20200309
parms$nStreams <- 1000

# Output and Progres report related
parms$outDir <- "../output/"
parms$start_time <- format(Sys.time(), "%Y%m%d_%H%M")
parms$report_file_name <- paste0("pooled_ZL2016-mc-", 
                                 parms$start_time, 
                                 ".txt")
parms$results_file_name <- paste0("pooled_ZL2016-mc-",
                                  parms$start_time,
                                  ".rds")
parms$description <- c("In each repetition, 1 dataset is created for each condition.
      Imputation methods are used on that condition-specific dataset.
      Results are therefore given per dataset in condition")

# Conditions --------------------------------------------------------------

p   <- c(200, 1000)   # number of variables
rho <- c(0)   # autoregressive structure
q   <- c(4)   # c(4, 20)  # active set (num of variables are true predictors of y)

conds <- as.matrix(expand.grid(p, rho, q))
