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
# library(rpart)     # for computing decision tree models
# library(treeClust) # for the rpart.predict.leaves(rp, newdata, type = "where") function
# library(dplyr)
library(randomForest)
library(CALIBERrfimpute) # package from Shah et al 2014

# Load all purpose functions
  source("./functions_allpurp.R")

# data --------------------------------------------------------------------
  # Create using datagen function
  source("./dataGen_test.R")
  set.seed(20191213)
  dt <- missDataGen(n=200, p=20)
    dt_c <- dt[[1]] # fully observed
    dt_i <- dt[[2]] # with missings
    dim(dt_i)



# Functions study ---------------------------------------------------------
  # Prep data
    K <- ncol(dt_i)-sum(tail(mice::md.pattern(dt_i),1) == 0) # number of variables needing imputation
    K_names <- names(which(colSums(apply(dt_i, 2, is.na)) != 0)) # select the names of these k variables
    Z <- dt_i    # dataset with missing values
    p <- ncol(Z) # number of variables [INDEX with j]
    n <- nrow(Z) # number of observations
    colnames(Z) <- paste0(rep("z_", p), seq(1, p))
    r <- colSums(apply(Z, 2, function (x) {!is.na(x)} )) # vector with number of observed values in a j variable

    # Define variables with missings
    l <- ncol(Z)-sum(tail(mice::md.pattern(Z),1) == 0) # number of variables needing imputation
    l_names <- names(which(colSums(apply(Z, 2, is.na)) != 0)) # select the names of these k variables
    
    # Define a dataset to recieve initialized missing values
    Z0 <- Z 
    
    # Define names of variables w/ missing values (general and by measurment scale)
    vartypes <- rbind(lapply(lapply(Z0[, names(Z0) %in% l_names], class), paste, collapse = " "))
    contVars <- colnames(vartypes)[vartypes == "numeric"] # names of continuous variables for selection
    factVars <- colnames(vartypes)[vartypes == "factor"] # includes factors and ordered factors
    ordeVars <- colnames(vartypes)[vartypes == "ordered factor"] # includes factors and ordered factors
    
    # Make oredered factors as numeric
    Z0[, ordeVars] <- as.numeric(Z0[, ordeVars])
    
    # Impute sample means for continuous variables and odered factors
    s_means <- apply(Z0[, c(contVars, ordeVars)], 2, mean, na.rm = TRUE) # sample means
    for (j in 1:length(c(contVars, ordeVars))) {
      Z0 <- Z0 %>% mutate_at(vars(c(contVars, ordeVars)[j]),
                             ~replace(., is.na(.), s_means[j])
      )
    }
    # Impute most common level for unordered factors
    for (j in 1:length(factVars)) {
      x <- addNA(Z0[, factVars[j]])
      m_commo <- names(which.max(table(x)))
      levels(x) <- c(levels(Z0[, factVars[j]]), 99)
      x[x == 99] <- m_commo
      x <- droplevels(x)
      Z0[, factVars[j]] <- x
    }
    
    Zm <- Z0
    
    
#### > mice.impute.rfcont() ####
    # Which variable to impute
    k <- 4 # for imputing one variable ()
    
    # inputs
    y <- Z[, l_names[k]]
    ry <- !is.na(y) # response variable T = observed, F = NA
    x <- Zm[, -which(colnames(Zm) == l_names[k])]
    ntree_cont    = 7    # number of trees
    nodesize_cont = NULL # minimum size of nodes
    maxnodes_cont = NULL # maximum number of nodes
    ntree = NULL         # 
    
    # function
    #x <- as.matrix(x) # mixed data does not become matrix easily
    bootsample <- sample(x=sum(ry), replace = TRUE) #sample takes place from 1:x
    
    # Random forest parameters
    if (is.null(ntree_cont) & !is.null(ntree)) {
      ntree_cont <- ntree
    }
    if (is.null(ntree_cont)) {
      if (is.null(getOption("CALIBERrfimpute_ntree_cont"))) {
        ntree_cont <- 10
      }
      else {
        ntree_cont <- getOption("CALIBERrfimpute_ntree_cont")
      }
    }
    if (is.null(nodesize_cont)) {
      if (is.null(getOption("CALIBERrfimpute_nodesize_cont"))) {
        nodesize_cont <- 5
      }
      else {
        nodesize_cont <- getOption("CALIBERrfimpute_nodesize_cont")
      }
    }
    if (is.null(maxnodes_cont)) {
      maxnodes_cont <- getOption("CALIBERrfimpute_maxnodes_cont")
    }
    if (ntree_cont > 1) { # SELECTION of BOOTSTRAPPED OBSERVATIONS
      yobs <- y[ry][bootsample]
      xobs <- x[ry, , drop = FALSE][bootsample, , drop = FALSE]
    } else {
      yobs <- y[ry]
      xobs <- x[ry, , drop = FALSE]
    }
    
    # Actual imputation
    xmiss <- x[!ry, , drop = FALSE]
    rf <- randomForest(xobs, yobs, ntree = ntree_cont, nodesize = nodesize_cont, 
                       maxnodes = maxnodes_cont)
    
    yhat <- predict(rf, xmiss,
                    predict.all=FALSE) # as means of the conditonal distributions
                                      # the values predicted using all of the random forest
                                      # trees are used (aggregate). If predict.all was FALSE
                                      # then, both aggregate and individual predictions (done
                                      # by each of the trees) would be reported.
    
    sqrtOOB <- sqrt(rf$mse[ntree_cont]) # for the last tree (just one of the 
                                        # ones we could have selected) the mean squared error of 
                                        # prediction made on the train data xobs
                                        # This does not make sense to me! Why the last tree? <- ????????????
    yimp <- rnorm(length(yhat), mean = yhat, sd = sqrtOOB)
    #return(yimp)
      
#### > mice.impute.rfcat() ####
    # Which variable to impute
    k <- 1 # for imputing one variable ()
    
    # inputs
    y <- Z[, l_names[k]]
    ry <- !is.na(y) # response variable T = observed, F = NA
    x <- Zm[, -which(colnames(Zm) == l_names[k])]
    ntree_cat = NULL
    nodesize_cat = NULL
    maxnodes_cat = NULL
    ntree = NULL
    
    if (is.logical(y)) {
      convertlogical <- TRUE
      y <- as.factor(y)
    } else {
      convertlogical <- FALSE
    }
    # Random forest parameters
    if (is.null(ntree_cat) & !is.null(ntree)) {
      ntree_cat <- ntree
    }
    if (is.null(ntree_cat)) {
      if (is.null(getOption("CALIBERrfimpute_ntree_cat"))) {
        ntree_cat <- 10
      }
      else {
        ntree_cat <- getOption("CALIBERrfimpute_ntree_cat")
      }
    }
    if (is.null(nodesize_cat)) {
      if (is.null(getOption("CALIBERrfimpute_nodesize_cat"))) {
        nodesize_cat <- 1
      }
      else {
        nodesize_cat <- getOption("CALIBERrfimpute_nodesize_cat")
      }
    }
    if (is.null(maxnodes_cat)) {
      maxnodes_cat <- getOption("CALIBERrfimpute_maxnodes_cat")
    }
    
    # Bootstrapping
    bootsample <- sample(sum(ry), replace = TRUE)
    yobs <- y[ry][bootsample]
    xobs <- x[ry, , drop = FALSE][bootsample, , drop = FALSE]
    xmiss <- x[!ry, , drop = FALSE]
    
    # Some stuff to fix the labeling (not imp)
    missinglevels <- (table(yobs) == 0)
    newlevels <- rep(NA_integer_, length(levels(y)))
    newlevels[!missinglevels] <- 1:sum(!missinglevels)
    labels <- levels(y)
    oldlevels <- 1:length(levels(y))
    changes <- !identical(newlevels, 1:length(levels(y)))
    if (changes) {
      temp <- data.frame(id_yobs = 1:length(yobs), fac = as.integer(yobs))
      lookup <- data.frame(fac = oldlevels, new = factor(newlevels))
      temp <- merge(temp, lookup, all.x = TRUE)
      yobs <- temp[order(temp$id_yobs), "new"]
    }
    
    # perform the imputations
    trees <- lapply(1:ntree_cat, function(x) {
      if (length(unique(yobs)) == 1) {
        yobs[1]
      }
      else {
        randomForest(xobs, yobs, ntree = 1, nodesize = nodesize_cat, 
                     maxnodes = maxnodes_cat)
      }
    })
    
    # This does not work beacuse you are keeping the data.frame structure instead
    # of the matrix form.
    yimp <- apply(xmiss, 1, function(x) {
      thetree <- trees[[sample(ntree_cat, 1)]]
      if ("randomForest" %in% class(thetree)) {
        predict(thetree, x)
      }
      else {
        thetree
      }
    })
    # This is what should happen inside for all the rows of xmiss
    thetree <- trees[[sample(ntree_cat, 1)]] # random tree selected and used for prediction
    predict(thetree, xmiss[1,])
    
    if (changes) {
      temp <- data.frame(id_yimp = 1:length(yimp), fac = as.integer(yimp))
      lookup <- data.frame(fac = newlevels, old = factor(oldlevels))
      levels(lookup$old) <- labels
      temp <- merge(temp, lookup, all.x = TRUE)
      yimp <- temp[order(temp$id_yimp), "old"]
    }
    if (convertlogical) {
      yimp <- as.logical(yimp == "TRUE")
    }
    return(yimp)
    
#### > mice.impute.rf() ####
    imp <- mice(dt_i, meth = "rf", ntree = 10)
    

# CART to adapt to RF with Shah -------------------------------------------

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
    
    K <- ncol(dt_i)-sum(tail(mice::md.pattern(dt_i),1) == 0) # number of variables needing imputation
    K_names <- names(which(colSums(apply(dt_i, 2, is.na)) != 0)) # select the names of these k variables
    
    # Initiliaze by stochastic regression imputation
    #   using only the variables with missings and the variable predicting
    #   the missing (cannot use all variables). Different aporaches could be
    #   followed here, and since I'm following the approach of Reiter and similar
    #   papers that are focused on creating sythetic datasets, there is not a
    #   specific guideline on what to do in this approach.
    print("Data initialization by MICE stochastic regression")
    imp <- mice::mice(dt_i[,1:p_imod], m = 1)
    dt_init <- cbind(mice::complete(imp)[, names(dt_i) %in% K_names], dt_i[, !names(dt_i) %in% K_names])
    
    dt_aug <- dt_init # this data will be the previous dataset at the beginning of the loop and the new data at the end
    
    ##### INITIALIZE #####
    
    # Match paper (see note up top) notation
    Z <- dt_i    # dataset with missing values
    p <- ncol(Z) # number of variables [INDEX with j]
    n <- nrow(Z) # number of observations
    colnames(Z) <- paste0(rep("z_", p), seq(1, p))
    r <- colSums(apply(Z, 2, function (x) {!is.na(x)} )) # vector with number of observed values in a j variable
    
    # For a jth variable say (example)
    j = 2
    z_j_obs   <- Z[, j][!is.na(Z[, j])]  # observed components of j-th variable
    z_j_mis   <- Z[, j][is.na(Z[, j])]   # missing components of j-th variable
    z_mj      <- Z[, -j]                 # collection of the p − 1 variables in Z except z_j (m for minus)
    z_mj_obs  <- z_mj[!is.na(Z[, j]), ]  # z_mj compponents corresponding to z_j_obs
    z_mj_mis  <- z_mj[is.na(Z[, j]), ]   # z_mj compponents corresponding to z_j_mis
    
    # Define variables with missings
    l <- ncol(Z)-sum(tail(mice::md.pattern(Z),1) == 0) # number of variables needing imputation
    l_names <- names(which(colSums(apply(Z, 2, is.na)) != 0)) # select the names of these k variables
    
    # Initiliaze by sample mean of variable with missing values
    print("Data initialization by sample mean")
    
    # Define a dataset to recieve initialized missing values
    Z0 <- Z 
    
    # Define names of variables w/ missing values (general and by measurment scale)
    vartypes <- rbind(lapply(lapply(Z0[, names(Z0) %in% l_names], class), paste, collapse = " "))
    contVars <- colnames(vartypes)[vartypes == "numeric"] # names of continuous variables for selection
    factVars <- colnames(vartypes)[vartypes == "factor"] # includes factors and ordered factors
    ordeVars <- colnames(vartypes)[vartypes == "ordered factor"] # includes factors and ordered factors
    
    # Make oredered factors as numeric
    Z0[, ordeVars] <- as.numeric(Z0[, ordeVars])
    
    # Impute sample means for continuous variables and odered factors
    s_means <- apply(Z0[, c(contVars, ordeVars)], 2, mean, na.rm = TRUE) # sample means
    for (j in 1:length(c(contVars, ordeVars))) {
      Z0 <- Z0 %>% mutate_at(vars(c(contVars, ordeVars)[j]),
                             ~replace(., is.na(.), s_means[j])
      )
    }
    # Impute most common level for unordered factors
    for (j in 1:length(factVars)) {
      x <- addNA(Z0[, factVars[j]])
      m_commo <- names(which.max(table(x)))
      levels(x) <- c(levels(Z0[, factVars[j]]), 99)
      x[x == 99] <- m_commo
      x <- droplevels(x)
      Z0[, factVars[j]] <- x
    }
    
    Zm <- Z0  # the dataset at Zm iteration 
              # when dataset hass been itialized, m=0, so Z0 = {z0_1, z0_2, ... , z0_l, z_l+1, ... , z_p}
              # each z_j of this data will be the m-1 "previous iteration" version at the beginning of 
              # the variable loop (for j in 1:l) and the current iteration data at the end
    
    ####
    
    for(i in 1:iters) {
      print(paste0("Iteration status: ", i, " of ", iters))
      for (k in 1:K) {
        print(paste0("Variable under imputation: ", K_names[k], " (", k, " of ", K, ")" ))
        # Prep Data
        indx_k <- names(dt_aug) %in% K_names[k]               # to index x_k
        ry <- !is.na(dt_i[, K_names[k]])                      # response variable T = observed, F = NA
        j <- 4
        Wm_j  <- Zm[,-j]
        zm_j  <- Zm[,j]
        Wm_mj <- Wm_j[is.na(Z[, j]), ]
        zm_mj <- zm_j[is.na(Z[, j])]
        
        ## Generate bootstrap sample
        indx_boSample <- sample(1:nrow(Wm_j), nrow(Wm_j),replace = TRUE)
        W.star.m_j    <- Wm_j[indx_boSample, ]
        z.star_j      <- zm_j[indx_boSample]
        

        
        rf_out <- randomForest(z.star_j ~ ., data=cbind(z.star_j, W.star.m_j))
        rf_out$predicted
        
        iris.rf <- randomForest(Species ~ ., data=iris, importance=TRUE, proximity=TRUE)
        
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