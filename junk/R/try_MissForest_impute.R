### Title:    imputeHD-comp impute w/ MissForest package
### Author:   Edoardo Costantini
### Created:  2019-DEC-13
### Modified: 2019-DEC-13
### Notes:    reference paper is Stekhoven Buhlmann 2012

# load packages
library(tidyverse)
library(dplyr)
library(mice)      
library(doParallel)
library(missForest)

# Prep data ---------------------------------------------------------------
# Load all purpose functions
  source("./functions_allpurp.R")
# Create using datagen function
  source("./dataGen_test.R")
    set.seed(20191120)
  dt <- missDataGen(n=500, p=1e3)
  dt_c <- dt[[1]] # fully observed
  dt_i <- dt[[2]] # with missings
    dim(dt_i)
    mice::md.pattern(dt_i)

# Impute w/ missForest ----------------------------------------------------

  set.seed(20191213)
  cl <- makePSOCKcluster(detectCores()-1) # to parallelize
  registerDoParallel(cl)
  system.time(
    
    #### IMPUTATION ####
    missForest_out <- missForest(dt_i,
                                 parallelize = "variable")
    ####################
    
  )
  stopCluster(cl)
  missForest_out$ximp