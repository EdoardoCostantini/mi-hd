### Title:    Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-05-19

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
  library(PcAux)
  library(caret) # for Elastic cv functions
  library(missForest) # for missForest imputation
  library(truncnorm) # for truncnorm sampling in probit models
  library(lavaan)

# Support Functions -------------------------------------------------------

  source("./functions.R")
  source("./functions_PCA_impute.R")
  source("./fun_DURR_impute.R")
  source("./fun_IURR_impute.R")
  source("./fun_blasso_impute.R")
  source("./fun_PCA_impute.R")
  source("./fun_MICE_CART_impute.R")
  source("./fun_MICE_RF_impute.R")
  source("./fun_MICE_TR_impute.R")
  source("./fun_RANF_impute.R")
  source("./fun_CART_impute.R")
  source("./fun_missFor_impute.R")
  source("./functions_genDt.R")
  source("./subroutines.R")

# Fixed Parameters --------------------------------------------------------

  parms <- list()

# Itereations, repetitions, etc
  parms$dt_rep     <- 10  # replications for averaging results (200 goal)
  parms$chains     <- 1 # number of parallel chains for convergence check
  parms$iters      <- 2 #50
  parms$burnin_imp <- 0 # 20 # how many imputation iterations should be discarded
  parms$ndt        <- 2 # 10 # number of imputed datasets to pool esitmaes from (10)
  parms$thin       <- (parms$iters - parms$burnin_imp)/parms$ndt
    # every how many iterations should you keep the imputation for a dataset
    # Example: of 20 iterations, I burn the first 10 I need for convergence
    # then, I keep 1 set of imputations every "thin"
    # number of iterations for the MI procedure (20)
  parms$pos_dt  <- (parms$burnin_imp+1):parms$iters # candidate datasets (after convergence)
  parms$keep_dt <- parms$pos_dt[seq(1, length(parms$pos_dt), parms$thin)] # keep 1 dataset every thin
  
  parms$rfntree <- 10  # default value as no indication is given in paper

  # For mice-like algorithms
  parms$mice_iters <- 2 #  20
  parms$mice_ndt   <- parms$ndt # 10 # number of imputed datasets to pool esitmaes from (10)
  
# Data gen ----------------------------------------------------------------

  parms$n       <- 200 # number of cases for data generation
  parms$cat     <- FALSE # categorical data or not? TRUE=yes, FALSE=no

# "True" values
  parms$item_mean <- 5
  parms$item_var  <- 5
  
# Variable blocks
  parms$blck1 <- 1:5 # higly correlated
  parms$blck2 <- (tail(parms$blck1,1)+1):10 # correlated

  parms$blck1_r <- .6 # correlation for highly correlated variables
  parms$blck2_r <- .3 # correlation for correlated variables

# Z gen (fully observed covariates)
  parms$Z_o_mu  <- 0 # mean for gen of fully observed variables

# z_m gen (covariates that will have missingness)
  parms$z_m_id  <- c((parms$blck1[1]):3, 
                     (parms$blck2[1]):8)
  parms$zm_n <- length(parms$z_m_id)
  parms$z_m_var <- 4 # error variance for z_m generation
  parms$S_all   <- list(q1 = (c(parms$blck1, parms$blck2)), # variables in block 1 and 2
                        q2 = (c(4:13, 50:59)))
  parms$stnr    <- c(1, sqrt(.2)) # signal to noise ratio for alpha coefs

# y gen / imporntant predictors
  parms$formula <- paste0("z", parms$z_m_id, collapse = ", ")
  # In the current set up there is no y, no formula, there are only
  # vairbales with missing values that are correlated with each
  # other. SO I use this object to identify the variables that
  # are to be imputed
  # parms$formula <- c("y ~ z1 + z2 + z3 + z4 + z5")
  # parms$k       <- 5 # number of predictors for analysis model

# Response Model (rm)
  parms$rm_b <- c(-3, -3, 2, -1)

  # parms$rm_x <- t(sapply(parms$z_m_id, function(x)
  #   c(base::sample(parms$blck1[which(parms$blck1 != x)], 3),
  #     base::sample(parms$blck2[which(parms$blck2 != x)], 1))
  # ))
  parms$rm_x <- matrix(c(2, 3, 5, 10,
                         5, 3, 1, 7,
                         1, 2, 5, 6,
                         3, 1, 2, 9,
                         2, 5, 1, 10,
                         5, 2, 1, 6), 
                       ncol = 4, nrow = 6,
                       byrow = TRUE)
  
  
# Models ------------------------------------------------------------------
  # source("./gen_lavaan_model.R") # generate txt file for lavaan model
  # # parms$lav_model <- read.table("../txt/lavaan_model_sat.txt",
  # #                         as.is = TRUE)$V1
  parms$lav_model <- paste(readLines("../txt/lavaan_model_sat.txt"), collapse="\n")
  
  
# Generic
  parms$meth_sel <- data.frame(DURR_la = TRUE,
                               DURR_el = FALSE,
                               IURR_la = TRUE,
                               IURR_el = FALSE,
                               blasso  = TRUE,
                               MI_PCA  = TRUE,
                               MI_CART = TRUE,
                               MI_RF   = TRUE,
                               MI_T    = TRUE,
                               missFor = TRUE,
                               GS      = TRUE,
                               CC      = TRUE
                               )
  parms$methods <- names(parms$meth_sel)[which(parms$meth_sel==TRUE)]
    # (GS, CC always last, alwyas present)

  parms$alphaCI <- .95 # confidence level for parameters CI
  parms$k_IURR  <- 0 # k value to bias coef sampling covariance matrix in IURR
                     # procedure to solve possible issues of singularity

# PCA method related
  parms$PCA_inter    <- FALSE # whether you want two way variables interactions
  parms$PCA_poly     <- FALSE # whether you want poly terms 
  parms$PCA_mincor   <- .3 # mincor for qucikpred for single imputation auxiliary vars
  parms$PCA_maxpw    <- 2L # polynomials order
  parms$PCA_pcthresh <- .5 # proportion of vairance for selection of PCs

# Replicability related
  parms$seed     <- 20200512 #20200309
  parms$nStreams <- 1000

# Output and Progres report related
  parms$outDir <- "../output/"
  parms$start_time <- format(Sys.time(), "%Y%m%d_%H%M")
  parms$report_file_name <- paste0("sim_res_", 
                                   parms$start_time, 
                                   ".txt")
  parms$results_file_name <- paste0("sim_res_", 
                                    parms$start_time,
                                    ".rds")
  parms$description <- c("In each repetition, 1 dataset is created for each condition.
        Imputation methods are used on that condition-specific dataset.
        Results are therefore given per dataset in condition")

# Conditions --------------------------------------------------------------

  p   <- 50 # c(50, 500) # number of variables
  rho <- c(0.5) # autoregressive structure
  q   <- c(1) # c(4, 20)  # active set (num of variables are true predictors of y)
  latent <- c(FALSE, TRUE)[1]
  pm <- c(.1, .3)
  
  conds <- expand.grid(p, rho, q, latent, pm)
    colnames(conds) <- c("p", "rho", "q", "latent", "pm")
