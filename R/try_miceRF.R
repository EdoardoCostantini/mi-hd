### Title:    imputeHD-comp impute with sequential CART
### Author:   Edoardo Costantini
### Created:  2019-DIC-16
### Modified: 2019-DIC-16
### Notes:    Main implementation is based on Shah et al 2014. mice package has a default version
###           based on Doove et al 2014 that differes from the Shah et al 2014 paper (see methods 
###           notes for details)

# load packages
library(tidyverse)
library(randomForest)
library(mice)
library(miceImpHDv) # modified version of mice with different Random Forest Algorithm

# Load all purpose functions
source("./functions_allpurp.R")

# data --------------------------------------------------------------------
# Create using datagen function
  source("./dataGen_test.R")
  set.seed(20191213)
  dt <- missDataGen(n=200, p=20)
  dt_c <- dt[[1]] # fully observed
  dt_i <- dt[[2]] # with missings

#### > mice.impute.rf() ####
  imp <- mice::mice(dt_i, meth = "rf", ntree = 10)

#### > mice.impute.rf() my version ####
  dt_i[missing_type(dt_i)$ordeVars] <- apply(dt_i[missing_type(dt_i)$ordeVars], 2, as.numeric)
    # does not support ordered factors as such, needs to be transfomred into numeric variables
  imp_rf_HDv <- miceImpHDv::mice(dt_i, meth = "rf", ntree = 10)
  
  # Compare results
  rf1 <- mice::complete(imp, "all")[[1]][is.na(dt_i[,1]), 1]
  rf2 <- mice::complete(imp, "all")[[2]][is.na(dt_i[,1]), 1]
  rf3 <- mice::complete(imp_rf_HDv, "all")[[1]][is.na(dt_i[,1]), 1]
  table(rf1==rf2)
  table(rf1==rf3)
  table(rf2==rf3)