### Title:    imputeHD-comp impute with bart
### Author:   Edoardo Costantini
### Created:  2019-NOV-14
### Modified: 2019-NOV-14
### Notes:    testing bart to impute a dataset

# load packages
library(tidyverse)
library(mice)
library(rpart)     # for computing decision tree models
library(dplyr)

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
#   the missing (cannot use all variables). Different aporaches could be
#   followed here, and since I'm following the aporach of Reiter and similar
#   papers that are focused on creating sythetic datasets, there is not a
#   specific guideline on what to do in this approach.
  
  imp <- mice::mice(dt_i[,1:7], method = "norm.predict", m = 1)
  dt_inizialised <- cbind(mice::complete(imp)[, names(dt_i) %in% K_names], dt_i[, !names(dt_i) %in% K_names])

# 
  for (k in 1:K) {
    k <- 1
    Xno_k <- dt_inizialised[, !names(dt_inizialised) %in% K_names[k]] # I want to use Xobs AND Xmiss(t-1) for building the trees
    x_k <- dt_inizialised[, names(dt_inizialised) %in% K_names[k]]  
      if(length(unique(x_k)) != 2) treetype <- "anova" else treetype <- "class" # anova is the option for regression tree
    data4tree <- data.frame(x_k=x_k, Xno_k)
    
    # Grow tree
    tree <- rpart(x_k ~., 
                       data = data4tree, 
                       method = treetype)
    plot(tree)
    tree$frame
    text(tree)
    
    # Identify who to predict
    Xobs <- dt_inizialised[is.na(dt_i[ , names(dt_inizialised) %in% K_names[k]]), colnames(X)]  # observed Xno_k for the missing on x_k
    
    # Identify the predictions
    leaf_id <- tree$where # terminal nodes
    leafs <- data.frame(x_k = dt_inizialised[, names(dt_inizialised) %in% K_names[k]], 
                        lf = tree$where)
    #arrange(leafs, lf)
    
    sample( leafs[leaf_id==6, 1], 1) # make this a bayesian bootstrapping sampling
    
  }