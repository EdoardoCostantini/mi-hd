### Title:    imputeHD-comp impute with bart
### Author:   Edoardo Costantini
### Created:  2019-NOV-13
### Modified: 2019-NOV-13
### Notes:    testing bart to impute a dataset

# load packages
library(BayesTree)
library(tidyverse)

# Prep data ---------------------------------------------------------------

  # Load
  dt <- readRDS("../data/data_tryout.rds") # generate with dataGen_test.R if not available
  dt_c <- dt[[1]] # fully observed
    mice::md.pattern(dt_c)
    dim(dt_c)
  dt_i <- dt[[2]] # with missings
    mice::md.pattern(dt_i)
    dim(dt_i)
    
  # Define variables with missings
  K <- ncol(dt_i)-sum(tail(mice::md.pattern(dt_i),1) == 0) # number of variables needing imputation
  K_names <- names(which(colSums(apply(dt_i, 2, is.na)) != 0)) # select the names of these k variables
  
# Imputation --------------------------------------------------------------
# Initialize by regression imputation
#   using only the variables with missings and the variable predicting
#   the missing (cannot use all variables).
  
  imp <- mice::mice(dt_i[,1:7], method = "norm.predict", m = 1)
  dt_imputed <- cbind(mice::complete(imp)[, names(dt_i) %in% K_names], dt_i[, !names(dt_i) %in% K_names])
    dim(dt_imputed)

# Initialize by sample mean (as suggested in Xu et al 2016 example 1)
  sample_means <- apply(dt_i[, names(dt_imputed) %in% K_names], 2, mean, na.rm = TRUE)
  dt_imputed <- dt_i
  for (k in 1:K) {
    dt_imputed <- dt_imputed %>% dplyr::mutate_at(vars(K_names[k]), ~replace(., is.na(.), sample_means[k]))
  }

# Imputations
  X0 <- dt_c[, !names(dt_c) %in% K_names] # set of covariates fully observed
  X1p <- dt_c[, names(dt_c) %in% K_names] # set of covariates with potentially missing values (I also include y)
  for (k in 1:K) {
    k <- 1
    X <- dt_imputed[, !names(dt_imputed) %in% K_names[k]] # I want to use Xobs AND Xmiss(t-1) for building the trees
    y <- dt_imputed[, names(dt_imputed) %in% K_names[k]]  # I want to use Xobs AND Xmiss(t-1) for building the trees
    Xobs <- dt_imputed[is.na(dt_i[ , names(dt_imputed) %in% K_names[k]]), colnames(X)]  # observed Xs for the missing ys
    bt<- bart(X, y, Xmiss, # using previously imputed dataset
              # general setup
              ntree = 10,
              ndpost = 1e3,
              keepevery = 1e3/2,
              # prior specifications
              # prior for error variance prior guess and degrees of freedom
              sigest = 1,
              sigdf = 3,      # nu
              sigquant = .90, # q
              # prior for the mean
              k = 2,          # varaince scale parameter of the prior for mean
              # prior for the trees
              power = 2, base = .95)
    ymiss_t <- t(bt$yhat.test)[,1] # I only keep the draw in the middle of the iterations (the first one)
  }
  
  
  