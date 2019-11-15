### Title:    imputeHD-comp impute with bart
### Author:   Edoardo Costantini
### Created:  2019-NOV-14
### Modified: 2019-NOV-15
### Notes:    Main implementation is based on Burgette Reiter 2010 Sequential Regression Trees
###           and Drechsler Reiter 2011
###           - Homogeneity criterion: Gini index for categorical, deviance for numerical
###             (look into source code to see what rpart uses for these)
###           - initialization stochastic regression imputation through mice package
###             following Burgette Reiter 2010 (draws from predictive distribution conditional on observed data)
###           - impute by sampling from a bayesian bootstrap sample of a leaft node (suggested by Drechsler Reiter 2011, 
###             implemented thanks to Rubin 1998 The Bayesian Bootstrap)

# load packages
library(tidyverse)
library(mice)
library(rpart)     # for computing decision tree models
library(treeClust) # for the rpart.predict.leaves(rp, newdata, type = "where") function
library(dplyr)
library(bayesboot)

# data --------------------------------------------------------------------

# load data
  dt <- readRDS("../data/data_tryout.rds") # generate with dataGen_test.R if not available
  dt_c <- dt[[1]] # fully observed
  mice::md.pattern(dt_c)
  dim(dt_c)
  dt_i <- dt[[2]] # with missings
  mice::md.pattern(dt_i)
  dim(dt_i)
  dt_i[1:10, 1:6]

# Define variables with missings
  K <- ncol(dt_i)-sum(tail(mice::md.pattern(dt_i),1) == 0) # number of variables needing imputation
  K_names <- names(which(colSums(apply(dt_i, 2, is.na)) != 0)) # select the names of these k variables

# Imputation w/ rpart -----------------------------------------------------

  output <- CARTimpute(dt_i, 10, 7)
    
# Functions ---------------------------------------------------------------

# Data prep
  arrange_dt_x_k <- function(dt, var_name){
  ## Description
    # Input: 
    # - dt: a complete dataset (initialized or augmented) containing a variable to be imputed and all imputation model variables
    # - var_name: character vector containing only the name of variable under imputation
    # Output: a dataset ready to be passed in a tree growing function with the under imputation variable in first column
  ## Trial inputs
    # dt <- mtcars
    # var_name <- names(mtcars)[3]
    indx_k <- names(dt) %in% var_name         # indexing obj: which variable is under imputation
    indx_no_k <- !names(dt) %in% var_name     # indexing obj: all vairables EXCEPT the one under imputation 
    Xno_k <- dt[, indx_no_k]                  # Use augmented Xno_k 
    x_k <- dt[, indx_k]                       # Select variable under imputation
    output <- data.frame(x_k = x_k, Xno_k)    # Combine the two for growing the tree
                                              # (doubt: should I grow the tree only on the observed values?)
    return(output)
  }
# Tree related  
  tree_type <- function(x){
  ## Description
    # Input: 
    # - x: a measured variable
    # Output: "class" / "anova" (this defines whether we need to grow a classification or a regression tree, respectively)
  ## Trial inputs
    # x <- mtcars$mpg # continuous example
    # x <- mtcars$am  # dichotomous example
    if(length(unique(x)) == 2) {type <- "class"} else {type <- "anova"}
    return(type)
  }
  
# Bayesian Bootstrap
  bbootstrap <- function(x, size) {
    # based on Rubin 1998
  ## Description
    # Input: 
    # - x: a variable of any type
    # - size: number of obseervations desired in the new bayesian bootstrap sample
    # Output: a bayesian bootsrap sample of size size of x
  ## Trial inputs
    # x <- rbinom(n, 1, .5)
    # size <- 10
    u <- sort(c(runif(length(x)-1, 0, 1))) # n-1 uniform draws
    g <- numeric(0)
    for (i in 1:(length(x))) {
      if(length(u[i-1]) == 0) u_prev <- 0 else u_prev <- u[i-1]
      g[i] <- u[i]-(u_prev)
      if(i == length(x)) {
        u[i] <- 1
        g[i] <- 1-(u[i-1])
      }
      #print(cbind(u[i], u_prev, g[i]) ) # check that it works
    }
    bbsample <- sample(x, 
                       size = size, 
                       replace = TRUE, 
                       prob = g)
    return(bbsample)
  }
  
# Single Imputation function
  CARTimpute <- function(dt_i, iter, p_imod) {
  ## Description
    # Input: 
    # - dt_i: a dataset with missing values
    # - iter: how many iterations one wants for the sequential cart algorithm
    # - p_imod: number of relevant predictors for imputation model (needed for stochastic regression for initial values)
    # Output: one imputed dataset
  ## Trial inputs
    # dt_i <- airquality #dt[[2]] # with missings mice::md.pattern(airquality)
    # iter <- 10
    # p_imod <- ncol(dt_i) # for low dimensional data this is enough 
    
    K <- sum(is.na(colSums(dt_i, na.rm = FALSE))) # number of variables needing imputation
    K_names <- names(which(colSums(apply(dt_i, 2, is.na)) != 0)) # select the names of these k variables
    
    # Initiliaze by stochastic regression imputation
    #   using only the variables with missings and the variable predicting
    #   the missing (cannot use all variables). Different aporaches could be
    #   followed here, and since I'm following the approach of Reiter and similar
    #   papers that are focused on creating sythetic datasets, there is not a
    #   specific guideline on what to do in this approach.
    imp <- mice::mice(dt_i[,1:p_imod], method = "norm.nob", m = 1)
    dt_inizialised <- cbind(mice::complete(imp)[, names(dt_i) %in% K_names], dt_i[, !names(dt_i) %in% K_names])
    
    dt_prev <- dt_inizialised # this data will be the previous dataset at the beginning of the loop and the new data at the end
    
    for(iter in 1:10) {
      for (k in 1:K) {
        # Prep Data
        indx_k <- names(dt_prev) %in% K_names[k]               # to index x_k
        data4tree <- arrange_dt_x_k(dt_prev, K_names[k])
        
        # Define tree type (classification or regression) based on nature of variable under imputation
        tt <- tree_type(data4tree$x_k)
        
        # Grow tree
        tree <- rpart(x_k ~., 
                      data = data4tree, 
                      method = tt)

        # Identify terminal nodes (l leafs) predictions (costum)
        leafs <- data.frame(x_k = dt_prev[, indx_k],
                            lf = tree$where)        # terminal nodes
        leaf_pred <- numeric(0)
        for (l in 1:length(unique(tree$where))) {
          x <- leafs[leafs$lf == unique(tree$where)[l] ,1] # for each terminal node, exctract the values of x_k in there
          bbx <- bbootstrap(x, length(x))                  # take a bootstrap sample of these values
          leaf_pred[l] <- sample(bbx, 1)                   # and sample 1 value from this conditional distribution
                                                           # (doubt: is a simple random sample with every value having equal prob enough?)
        }
        imputation_frame <- data.frame(leaf_id = unique(tree$where), # what value should be imputed if you fall in a node 
                                       leaf_pred)
        
        # Identify Xno_k values corresponding to the missing x_k
        Xobs <- dt_prev[is.na(dt_i[ , indx_k]), !indx_k]  # observed Xno_k for the missing on x_k
        
        # Impute values
        missings_frame <- rpart.predict.leaves(tree, Xobs, type = "where") # locate obs w/ missing x_K in a tree based on Xobs
        for (obs in 1:length(missings_frame)) {
          dt_prev[names(missings_frame)[obs], indx_k] <- imputation_frame[imputation_frame$leaf_id == missings_frame[obs], 2]
          #imputation_frame %>% filter(leaf_id == missings_frame[obs]) %>% select(2)
        }
      }
      #print(dt_prev %>% filter(is.na(dt_i$yinc)) %>% select(1:7))
    }
    return(dt_prev)
  }
  
# Imputation w/ MICE cart implementation ----------------------------------

  