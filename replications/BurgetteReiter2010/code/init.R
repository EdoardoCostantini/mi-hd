### Title:    Initialize Environment and Parameters for Imputed DV Simulation
### Author:   Anonymized for peer review
### Created:  2020-02-13
### Modified: 2020-02-15

## Optional: to perform checks

# library(mvtnorm)
# library(mice)   
# source("./functions.R") # applyLib()
# source("./subroutines.R")    
# source("./simMissingness.R") 

# packages ----------------------------------------------------------------

library(rlecuyer) # for set seeds

# Simulation set up -------------------------------------------------------

## Define the fixed simulation parameters:

parms <- list()

parms$n <- 1e3      # sample size of the datasets to be generated
parms$p <- 10       # count of predictors X (excludes y)
parms$b_true <- c(0, 0, .5,.5,.5,.5,.5,1,1) # true values for B
parms$b <- 9        # count of coefficient to estimate (includes intercept)
parms$n_meth <- 2   # count of imputation methods compared
parms$dt_rep <- 1e3 # count of datasets to average over
parms$chains <- 10  # count of times a dataset should be imputed (m)
parms$iters <- 10   # mice iterations in 1 chain
parms$seed <- 20200214
parms$nStreams <- 1000

