### Title:    Replication Zhao Long 2016 - initialization script
### Author:   Edoardo Costantini
### Created:  2020-02-20
### Modified: 2020-02-20

# Packages ----------------------------------------------------------------

library(tidyverse)
library(CVTuningCov) # for AR1 function
library(mvtnorm)
library(monomvn)
library(glmnet)
library(rlecuyer) # for seed in parallel
library(blasso) # for blasso function according to Hans2010
library(parallel)

# Support Functions -------------------------------------------------------

source("./functions.R")
source("./fun_DURR_impute.R")
source("./fun_IURR_impute.R")
source("./fun_blasso_impute.R")
source("./subroutines.R")

# Fixed Parameters --------------------------------------------------------

parms <- list()

# Itereations, repetitions, etc
parms$dt_rep  <- 500 # replications for averaging results (500 goal)
parms$chains  <- 10  # number of imputed datasets to pool esitmaes from
parms$iters   <- 1   # number of iterations for the MI procedure
parms$iter_bl <- 1e3 # blasso samples
parms$burn_bl <- 5 * parms$iter_bl # blasso burnin

# Data Gen
parms$n       <- 100 # number of cases for data generation

# Z gen (fully observed covariates)
parms$Z_o_mu  <- 0 # mean for gen of fully observed variables

# z_m gen (covariates that will have missingness)
parms$z_m_id  <- c(1, 2, 3) # variable ids will have missing values
parms$z_m_var <- 4   # error variance for z_m generation
parms$S_all   <- list(q4=(c(4, 5, 50, 51)-1),
                      q20=(c(4:13, 50:59)-1))
parms$stnr    <- c(1, sqrt(.2)) # signal to noise ratio for alpha coefs

# y gen
parms$y_var   <- 6   # error variance for y generation
parms$b       <- 1   # betas and intercept in linear model to gen y
parms$formula <- "y ~ z1 + z2 + z3 + z4 + z5" # analysis model
parms$k       <- 5 # number of predictors for analysis model

# impose missingness
parms$b_miss_model <- c(-1, -1, 2, -1)  # betas for generation of missingness
parms$detlamod <- list(c("z4", "z5", "y"), # predictors logit model miss impose
                       c("z4", "z51", "y"),
                       c("z50", "z51", "y"))

# Generic
parms$alphaCI <- .95 # confidence level for parameters CI
parms$k_IURR <- 0 # k value to bias coef sampling covariance matrix in IURR
                  # procedure to solve possible issues of singularity

parms$methods <- c("DURR", "IURR", "blasso", "MICE-RF", "MI_T", "CC") 
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
q   <- c(4, 20)   # c(4, 20)  # active set (num of variables are true predictors of y)

conds <- as.matrix(expand.grid(p, rho, q))
