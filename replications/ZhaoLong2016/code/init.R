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

# Support Functions -------------------------------------------------------

source("./functions.R")
source("./fun_DURR_impute.R")
source("./fun_IURR_impute.R")
source("./fun_BLasso_impute.R")
source("./subroutines.R")

# Fixed Parameters --------------------------------------------------------

parms <- list()
parms$dt_rep  <- 50 # 500 on original
parms$chains  <- 10
parms$iters   <- 1
parms$n       <- 100
parms$y_var   <- 3
parms$b       <- 1 # beta in linear model to create y (same for all x and intercept)
parms$b_miss_model <- c(-1, -.1, 2, -2)
parms$stnr    <- c(1, sqrt(.2), sqrt(.08)) # signal to noise ratio
parms$S_all   <- list(q4=(c(2,3,50,51)-1),
                      q20=(c(2:11, 50:59)-1),
                      q50=(c(2:11, 50:59, 70:79, 90:99, 110:119)-1))
parms$formula <- "y ~ z1 + z2 + z3"
parms$alphaCI <- .95 # confidence level for parameters CI
parms$k       <- 3 # number of predictors for analysis model
parms$iter_bl <- 40 # 1e3
parms$burn_bl <- 1/10 * parms$iter_bl #

parms$k_IURR <- 0 # k value to bias coef sampling covariance matrix in IURR
                  # procedure to solve possible issues of singularity

parms$methods <- c("DURR", "IURR", "Blas", "MI_T", "MI_50", "CC") 
  # for naming objects, hence order is important (CC always last)

parms$seed     <- 20200309
parms$nStreams <- 1000



# Conditions --------------------------------------------------------------

p <- c(200, 1000)   # number of variables
rho <- c(0) # autoregressive structure
q <- c(4, 20)   # active set (num of variables are true predictors of y)

conds <- as.matrix(expand.grid(p, rho, q))
