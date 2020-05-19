### Title:    Replication Deng Et Al 2016 - initialization script
### Author:   Edoardo Costantini
### Created:  2020-05-05
### Notes: With this init file you can formulate the set up for ZL2016 replication
###        and test that set up with the multivariate generalization of imputation 
rm(list=ls())

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
source("./fun_MICE_RF_impute.R")
source("./fun_MICE_TR_impute.R")
source("./subroutines.R")

# Fixed Parameters --------------------------------------------------------

parms <- list()

# Itereations, repetitions, etc
parms$dt_rep  <- 500 # replications for averaging results (200 goal)
parms$chains  <- 5 # number of parallel chains for convergence check
parms$ndt     <- 10 # number of imputed datasets to pool esitmaes from (10)
parms$iters   <- 1
parms$burnin_imp <- 10 # how many imputation iterations should be discarded
parms$thin    <- 1 
  # every how many iterations should you keep the imputation for a dataset
  # Example: of 20 iterations, I burn the first 10 I need for convergence
  # then, I keep 1 set of imputations every "thin"
  # number of iterations for the MI procedure (20)
parms$pos_dt <- (parms$burnin_imp+1):parms$iters # candidate datasets (after convergence)
parms$keep_dt <- parms$pos_dt[seq(1, length(parms$pos_dt), parms$thin)] # keep 1 dataset every thin

parms$iter_bl <- 1  # blasso samples (everything else is in burnin)
parms$burn_bl <- 1e3 * parms$iter_bl # blasso burnin
parms$rfntree <- 10 # default value as no indication is given in paper

# Data Gen
parms$n       <- 100 # number of cases for data generation

# Z gen (fully observed covariates)
parms$Z_o_mu  <- 0 # mean for gen of fully observed variables

# z_m gen (covariates that will have missingness)
parms$z_m_id  <- 1
parms$z_m_var <- 1 # error variance for z_m generation
parms$S_all   <- list(q4=(c(2,3,50,51)),
                           q20=(c(2:11, 50:59)))
parms$stnr    <- c(1, sqrt(.2)) # signal to noise ratio for alpha coefs

# y gen
parms$y_var   <- c(3, 6)[paper] # error variance for y generation
parms$b       <- 1   # betas and intercept in linear model to gen y
parms$formula <- c("y ~ z1 + z2 + z3")   # analysis model
parms$k       <- 3 # number of predictors for analysis model

# impose missingness
parms$b_miss_model <- c(-1, -.1, 2, -2) # betas for generation of missingness 

# predictors logit model miss impose
parms$detlamod <- c("z2", "z3", "y")

# Generic
parms$alphaCI <- .95 # confidence level for parameters CI
parms$k_IURR <- 0 # k value to bias coef sampling covariance matrix in IURR
                  # procedure to solve possible issues of singularity

parms$methods <- c("DURR", "IURR", "blasso", "MICE-RF", "MI_T", "GS", "CC")
  # for naming objects, hence order is important (CC always last)
  # c("DURR", "IURR", "blasso", "MICE_RF", "MI_T", "GS", "CC") # complete correct order

# Replicability related
parms$seed     <- 20200512 #20200309
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

p   <- c(200)
rho <- c(0)     # autoregressive structure
q   <- c(4, 20) # active set (num of variables are true predictors of y)

conds <- as.matrix(expand.grid(p, rho, q))
  colnames(conds) <- c("p", "rho", "q")
