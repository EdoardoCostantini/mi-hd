### Title:    imputeHD-comp impute w/ PCA auxiliary variables
### Author:   Edoardo Costantini
### Created:  2019-NOV-26
### Modified: 2019-NOV-26
### Notes:    Main implementation is based on Howard et al 2015
###           - Inclusion of auxiliary vairables: for each variable with missing values, a dataset
###             with only the variable and the first n PC components are use (n should be dedided,
###             look into the paper again)
###           - Initialization: auxiliary variables are initialized with one non stocastich regression
###             run. The idea is that we are not worried about standard error at tgis level. We simply 
###             want a full dataset for the PCA computation.

# load packages
library(tidyverse)
library(BayesTree)
library(tgp)      # for Bayesian CART implementation bcart
library(rpart)    # for computing decision tree models
library(bartpkg1) # github version of the package
library(sbart)    # newest cran version available

# Prep data ---------------------------------------------------------------

# Create using datagen function
source("./dataGen_test.R")
  set.seed(20191120)
dt <- missDataGen(n=100, p=30)
  dt_c <- dt[[1]] # fully observed
  dt_i <- dt[[2]] # with missings
    dim(dt_i)
    mice::md.pattern(dt_i)

# Define variables with missings
K <- ncol(dt_i)-sum(tail(mice::md.pattern(dt_i),1) == 0) # number of variables needing imputation
K_names <- names(which(colSums(apply(dt_i, 2, is.na)) != 0)) # select the names of these k variables

# Imputation --------------------------------------------------------------

  # Data preparetion
  # PCs extraction
    pr_out <- prcomp(dt_i[, !names(dt_i) %in% K_names], scale = TRUE)
    
    pr_var <- pr_out$sdev**2
    p_pr_var <- pr_var/sum(pr_var) # proportion of variance explained by component
  
  # Inlcude auxiliary PCs in the imputation
    dt_i[,1:5]
    dt_i_PCaux
    
  # Initialize dataset
  print("Data initialization by MICE stochastic regression")
  imp <- mice::mice(dt_i[,1:ncol(dt_i)], method = "norm.nob", m = 1, maxit = 1)
  dt_init <- cbind(mice::complete(imp)[, names(dt_i) %in% K_names], dt_i[, !names(dt_i) %in% K_names])
  
  dt_aug <- dt_init # this data will be the previous dataset at the beginning of the loop and the new data at the end
  