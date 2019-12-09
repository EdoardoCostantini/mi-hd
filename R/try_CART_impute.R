### Title:    imputeHD-comp impute with sequential CART
### Author:   Edoardo Costantini
### Created:  2019-NOV-14
### Modified: 2019-NOV-26
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

# Load all purpose functions
  source("./functions_allpurp.R")

# data --------------------------------------------------------------------
# Create using datagen function
source("./dataGen_test.R")
set.seed(20191120)
dt <- missDataGen(n=1e4, p=500)
dt_c <- dt[[1]] # fully observed
dt_i <- dt[[2]] # with missings
dim(dt_i)
mice::md.pattern(dt_i)

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
    # x <- as.factor(mtcars$am)  # dichotomous example
    if(is.factor(x)) {type <- "class"} else {type <- "anova"}
    return(type)
  }
  
# Bayesian Bootstrap
  bbootstrap <- function(x) {
    # based on Rubin 1998
  ## Description
    # Input: 
    # - x: a variable of any type
    # Output: a bayesian bootsrap sample of size size of x
  ## Trial inputs
    # x <- rbinom(30, 1, .5)
    # x <- rownames(airquality)
    size <- length(x)
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
    #sum(g)
    bbsample <- sample(x, 
                       size = size, 
                       replace = TRUE, 
                       prob = g)
    return(bbsample)
  }
  
# Single Imputation function
  CARTimpute <- function(dt_i, iters, p_imod) {
  ## Description
    # Input: 
    # - dt_i: a dataset with missing values
    # - iter: how many iterations one wants for the sequential cart algorithm
    # - p_imod: number of relevant predictors for imputation model (needed for stochastic regression for initial values)
    # Output: one imputed dataset
  ## Trial inputs
    # dt_i <- airquality # with missings mice::md.pattern(airquality)
    # iter <- 10
    # p_imod <- ncol(dt_i) # for low dimensional data this is enough 
    # k <- 1 # just for running the inside of the loop
  ## Notes: right now this implementation differes from Doove et al 2014 in that: (a) initialization of 
    # missing values is done with conditional mean imputation instead of random draws from the observed
    # relevant values; (B) I implemented a bayesian bootstrap resampling of the donor set following 
    # Burgette Rieter 2010; (C) I need to look more into the specification of a CART (do I need to define
    # the value of "cp")
    
    K <- sum(is.na(colSums(dt_i, na.rm = FALSE))) # number of variables needing imputation
    K_names <- names(which(colSums(apply(dt_i, 2, is.na)) != 0)) # select the names of these k variables
    
    # Initiliaze by stochastic regression imputation
    #   using only the variables with missings and the variable predicting
    #   the missing (cannot use all variables). Different aporaches could be
    #   followed here, and since I'm following the approach of Reiter and similar
    #   papers that are focused on creating sythetic datasets, there is not a
    #   specific guideline on what to do in this approach.
    print("Data initialization by MICE stochastic regression")
    imp <- mice::mice(dt_i[,1:p_imod], method = "norm.nob", m = 1)
    dt_init <- cbind(mice::complete(imp)[, names(dt_i) %in% K_names], dt_i[, !names(dt_i) %in% K_names])
    
    dt_aug <- dt_init # this data will be the previous dataset at the beginning of the loop and the new data at the end
    
    for(i in 1:iters) {
      print(paste0("Iteration status: ", i, " of ", iters))
      for (k in 1:K) {
        print(paste0("Variable under imputation: ", K_names[k], " (", k, " of ", K, ")" ))
        # Prep Data
        indx_k <- names(dt_aug) %in% K_names[k]               # to index x_k
        ry <- !is.na(dt_i[, K_names[k]])                      # response variable T = observed, F = NA
        data4tree <- arrange_dt_x_k(dt_aug, K_names[k])[ry, ] # arrange and select Xobs (values of Xno_k for which X_k is observed)
        
        # Define tree type (classification or regression) based on nature of variable under imputation
        tt <- tree_type(data4tree$x_k)
        
        # Grow tree
        tree <- rpart(x_k ~., data = data4tree, method = tt)
        
        # Identify terminal nodes (l leafs) predictions (costum)
        leafs <- data.frame(x_k = data4tree[, indx_k],
                            lf = tree$where)        # terminal nodes
        # Identify donor pool
        donors <- lapply(unique(tree$where), function(s) data4tree$x_k[tree$where == s]) # Create donor pools (for each node, select the y values 
                                                                               # that are observed (or previously imputed) for all cases in that node
        # Identify prediction for donor pool
        bbdonors <- lapply(donors, bbootstrap) # Bootstrap sample those donors
        leaf_pred <- vapply(seq_along(bbdonors), function(s) sample(bbdonors[[s]], 1), FUN.VALUE = numeric(1))
          # Get a simple sample from each donor pool that will be the imputed value for an observation that falls in a specific group
          # (doubt: is a simple random sample with every value having equal prob enough?)
        imputation_frame <- data.frame(leaf_id = unique(tree$where), # what value should be imputed if you fall in a leaf
                                       leaf_pred)
        # Xmis = values on X_no_k for observations with missing values on x_k
        Xmis <- dt_aug[!ry, !indx_k]
        Xmis_loca <- rpart.predict.leaves(tree, Xmis, type = "where") # location in the tree of Xmiss
        
        # Impute values
        for (obs in 1:length(Xmis_loca)) {
          dt_aug[names(Xmis_loca)[obs], indx_k] <- imputation_frame[imputation_frame$leaf_id == Xmis_loca[obs], 2]
        }
      }
      #print(dt_aug %>% filter(is.na(dt_i$yinc)) %>% select(1:7))
    }
    return(dt_aug)
  }
  
# Imputation w/ MICE cart implementation ----------------------------------
# for function details, see: https://github.com/stefvanbuuren/mice/blob/master/R/mice.impute.cart.R
  
# Use
  imp <- mice(dt_i, meth = 'cart', minbucket = 4)
  imp_sets <- mice::complete(imp)
  
# Sections of mice code that are relevant
    y <- dt_i$yinc
      ry <- !is.na(y)
      wy <- !ry
    x <- dt_i[, 10:50]

    cp <- 1e-04
    minbucket <- 5
    
    xobs <- data.frame(x[ry, , drop = FALSE])
    xmis <- data.frame(x[wy, , drop = FALSE])
    yobs <- y[ry]

    fit <- rpart(yobs ~ ., data = cbind(yobs, xobs), method = "anova", 
                 control = rpart.control(minbucket = minbucket, cp = cp))
    leafnr <- floor(as.numeric(row.names(fit$frame[fit$where, ])))          # leaf location of the training units
    fit$frame$yval <- as.numeric(row.names(fit$frame))                      # combined with the following predict finds the location 
    nodes <- predict(object = fit, newdata = xmis)                          # of the new data in the tree
    donor <- lapply(nodes, function(s) yobs[leafnr == s])
    impute <- vapply(seq_along(donor), function(s) sample(donor[[s]], 1), numeric(1))

  # Mice does not do the bootstrap sampling and does not initialize with conditional means