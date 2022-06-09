### Title:    Replication Howard et al 2015 - initialization script
### Author:   Anonymized for peer review
### Created:  2020-05-18

rm(list=ls())

# Packages ----------------------------------------------------------------

packs <- list("parallel",
              "rlecuyer", # for parallelisation
              "mvtnorm",
              "mice",
              "lavaan") # for FIML saturated model

lapply(packs, require, character.only = TRUE)

# library(tidyverse)
# library(CVTuningCov) # for AR1 function
# library(mvtnorm)
# library(monomvn)
# library(glmnet)
# library(rlecuyer) # for seed in parallel
# library(blasso) # for blasso function according to Hans2010
# library(parallel)

# Support Functions -------------------------------------------------------

source("./functions.R")
source("./subroutines.R")
source("./lavan_models.R")
# source("./fun_DURR_impute.R")
# source("./fun_IURR_impute.R")
# source("./fun_blasso_impute.R")
# source("./fun_MICE_RF_impute.R")
# source("./fun_MICE_TR_impute.R")

# Fixed Parameters --------------------------------------------------------

parms <- list()

# Itereations, repetitions, etc
parms$dt_rep  <- 1e2 #5e3 # replications for averaging results (200 goal)

# Data Gen
parms$vars <- c("Y", "X", paste0("A", seq(1, 8)))
parms$n_vars <- length(parms$vars)
parms$n <- 1e5
parms$sd <- 1
parms$mean <- 0
parms$cor_fix <- .3 # fixed correlation XY and X and A-1-8

parms$cor_rng <- c(.3, .6) #seq(.1, .6, by = .1) # range of YX correlation values
parms$mis_rate <- .6 #seq(.1, .8, by = .1) # rate of missing data
parms$N <- c(75, 200, 800, 1e3) # c(1e5, seq(50, 1e3, by = 25)) # sa,èòe soze
parms$mis_mech <- c(l=1L, nl=2L)[1] # linear or nonlinear MAR mechanism
parms$methods <- c("NAV", "IAV", "PCA", "GSp", "GSs") # imputation methods
parms$out_ref <- c(XY = "Y~~X", mu_y = "Y~1", var_y = "Y~~Y") # parameters of interest

# Replicability related
parms$seed     <- 20200512 #20200309
parms$nStreams <- 1000

# Output and Progres report related
parms$outDir <- "../output/"
parms$start_time <- format(Sys.time(), "%Y%m%d_%H%M")
parms$report_file_name <- paste0("pooled_HEAT2015_",
                                 parms$start_time,
                                 ".txt")
parms$results_file_name <- paste0("pooled_HEAT2015_",
                                  parms$start_time,
                                  ".rds")
parms$description <- c("In each repetition, 1 dataset is created for each condition.
      Imputation methods are used on that condition-specific dataset.
      Results are therefore given per dataset in condition")

# Condition matrix --------------------------------------------------------

conds <- as.matrix(expand.grid(YA_rho = parms$cor_rng, 
                               miss_r = parms$mis_rate, 
                               N      = parms$N, 
                               misMch = parms$mis_mech))

conds <- conds[2, ]