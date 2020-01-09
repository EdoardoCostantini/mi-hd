### Title:    imputeHD-comp impute with bart
### Author:   Edoardo Costantini
### Created:  2019-NOV-13
### Modified: 2019-DEC-12
### Notes:    testing bart to impute a dataset

# load packages
library(tidyverse)
library(BayesTree)
library(tgp)      # for Bayesian CART implementation bcart
library(rpart)    # for computing decision tree models
library(bartpkg1) # github version of the package
library(sbart)    # newest cran version available
library(BART)

# Prep data ---------------------------------------------------------------
  source("./functions_allpurp.R")

  # Create using datagen function
  source("./dataGen_test.R")
  set.seed(20191120)
  dt <- missDataGen(n=25, p=10)
  dt_c <- dt[[1]] # fully observed
  dt_i <- dt[[2]] # with missings
    mice::md.pattern(dt_i)
  
  # Define variables with missings
  r <- !is.na(dt_i)
  missVarInfo <- missing_type(dt_i)
  l <- missVarInfo$l
  l_names <- missVarInfo$l_names
  
# Imputation --------------------------------------------------------------
# Initialize by regression imputation
#   using only the variables with missings and the variable predicting
#   the missing (cannot use all variables).
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

# Initialize by sample mean (as suggested in Xu et al 2016 example 1)
  print("Data initialization by sample mean")
  sample_means <- apply(dt_i[, names(dt_i) %in% K_names], 2, mean, na.rm = TRUE)
  dt_aug <- dt_i
  for (k in 1:K) {
    dt_aug <- dt_aug %>% dplyr::mutate_at(vars(K_names[k]), ~replace(., is.na(.), sample_means[k]))
  } 
  head(dt_aug[,1:7], 10)
  head(dt_i[,1:7], 10)
# Sequential BART
  iters <- 10
  fast_diagnostic <- matrix(rep(NA, iters*7), ncol = 7)
for(i in 1:iters) {
  #i <- 1
  print(paste0("Iteration status: ", i, " of ", iters))
  for (k in 1:K) {
    #k <- 1
    print(paste0("Variable under imputation: ", K_names[k], " (", k, " of ", K, ")" ))
    # Prep Data
    indx_k <- names(dt_aug) %in% K_names[k]               # to index x_k
    ry <- !is.na(dt_i[, K_names[k]])                      # response variable T = observed, F = NA
    data4tree <- arrange_dt_x_k(dt_aug, K_names[k])[ry, ] # arrange and select Xobs (values of Xno_k for which X_k is observed)
    Xobs <- data4tree[, -1]                               # after rearranging, the variable under imputation is always first
      dim(Xobs)
    yobs <- data4tree[, 1] #the x_k, variable under imputation
      length(yobs)
    Xmiss <- arrange_dt_x_k(dt_aug, K_names[k])[!ry, -1]     # observed and previously imputed / initialized values on all other variables
      dim(Xmiss)
    # Define tree type (classification or regression) based on nature of variable under imputation
    #tt <- tree_type(data4tree$x_k)
    
    # Sample tree from distribution of trees given augmented data (step 1 in algorithm p.593)
    set.seed(1234)
    bt <- bart(Xobs, yobs, Xmiss, # using previously imputed dataset
               # general setup
               ntree = 50,
               ndpost = 1e3,
               keepevery = 1,
               # prior specifications
               # prior for error variance prior guess and degrees of freedom
               sigest = 1,
               sigdf = 3,      # nu
               sigquant = .90, # q
               # prior for the mean
               k = 2,          # varaince scale parameter of the prior for mean
               # prior for the trees
               power = 2, base = .95,
               verbose = F)
    bt$yhat.test # for observation 1, column 1 has many ndpost/keeprecovery values.
                 # each value is a sum of the mu value found by dropping that case in
                 # the m trees that have been grown, + some error. Each value is then 
                 # one sampled from the conditional distribution of Y given X, for a 
                 # specific X.
    bt$yhat.test.mean
                 # The mean of what has been described for yhat.test is the estiamte 
                 # of the Expected value of Y given that X takes the specific values
                 # of one observation for which you want a prediction. Hence, this value
                 # is a reasonable prediction of Y at a particular X
    # Impute values
    # Sample one mising value at the time
    y_miss_k <- NULL
    i <- 1
    # Sample proposition
    y_miss_k[i] <- rnorm(1, 
                         tail(bt$yhat.test[, i], 1), # I'm using the last draw from all the ones obtained with BART
                                                     # I'm not sure this is the correct procedure but for now let's stick to it
                         tail(bt$sigma, 1))
    # Accept/reject proposition
    data4tree[i,1] <- y_miss_k[i]
    full_data <- data.frame(data4tree[,-1], 
                            focusK = data4tree[,1])
    BCART <- rpart(focusK ~., data = full_data) # this should be a Bayesian tree as in tgp::bcart
    BCART$
    head(Xobs[, 1:7])
    dnorm(Xobs[i, k+1], )
    dt_aug[row.names(Xmiss), indx_k] <- bt$yhat.test.mean
  }
   new_impute <- dt_aug %>% filter(is.na(dt_i$yinc)) %>% select(1:7)
   original <- dt_c %>% filter(is.na(dt_i$yinc)) %>% select(1:7)
   fast_diagnostic[i,] <- colSums(new_impute - original)
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
  

# using 'bartpkg1' package -----------------------------------------------

  # Create dataset
  # Create using datagen function
  source("./dataGen_test.R")
  set.seed(20191120)
  dt <- missDataGen(n=50, p=7)

  dt_i <- data.frame(dt$complete$y, dt$incomplete[,-1])
    names(dt_i) = names(dt[[1]])
    # sbart implementatio requires a fully observed y 
    # and a covariates matrix with missing values
  mice::md.pattern(dt_i)
  dim(dt_i)
  
  # Impute
  imputedList <- serBARTfunc(dt_i[,-1], dt_i$y, 
                             datatype=0, #continuous covariates
                             type = 0,   #continuous dependent
                             numskip = 199, burn = 1000, m = 200,
                             sigdf = 3, sigquant = 0.9, kfac = 2)
  str(imputedList)

# using 'sbart' package ---------------------------------------------------
  imputedList <- seqBART(dt_i[,-1], dt_i$y, 
                         datatype=0, #continuous covariates
                         type = 0,   #continuous dependent
                         numskip = 199, burn = 1000, m = 200,
                         sigdf = 3, sigquant = 0.9, kfac = 2)
  str(imputedList)
  


  
# using bart directly -----------------------------------------------------
  # Initialize Dataset
  Z0 <- init_dt_i(dt_i, missing_type(dt_i))
  
  # This takes place for one variable, for one itaration, for one imputated dataset
  j <- 1 # define variable to be imputed
  l_names[[j]]
  
  # define y, ry and wy
  y  <- Z0[, j]
  x  <- model.matrix(Z0[, j]~., Z0[,-j])

  ry <- r[, j] # TRUE = observed
  wy <- !ry    # TRUE = NA
  
  # select cases to train and predict
  xobs <- x[ry, , drop = FALSE]
  xmis <- x[wy, , drop = FALSE]
  yobs <- y[ry]
  if(tree_type(yobs)=="class") {
    originallvls <- levels(yobs)
    yobs <- as.numeric(yobs)
    yobs[yobs == 1] <- 0
    yobs[yobs == 2] <- 1
    tt <- "pbart"
  } else {tt <- "wbart"}
  
  # grow bart
  m <- 200
  ndpost <- 15
    rb = gbart(x.train = xobs,
               y.train = yobs,
               x.test  = xmis, 
               # general setup
               type = tt,
               ntree = m,
               ndpost = ndpost,
               nskip = 1e3, # iterations treated as burn in
               # prior specifications
               # prior for error variance prior guess and degrees of freedom
               sigest = 1,
               sigdf = 3,      # nu
               sigquant = .90, # q
               # prior for the mean
               k = 2,          # varaince scale parameter of the prior for mean
               # prior for the trees
               power = 2, base = .95)  # power and base parameter for the tree prior
    
    indx <- 5
    plot(density(rb$yhat.test[,indx]))
    plot(seq(1, ndpost), rb$yhat.test[,indx], type = "l")
    
    rb$sigma         # ndpost draws of sigma
    rb$yhat.train    # ndpost draws of p(f | data)
    rb$yhat.test     # ndpost draws of predictive distribution
    rb$yhat.test[,1] # ndpost draws from posterior predictive distribution for observation 1, 
                     # of variable under imputation. Each prediction is based on a different
                     # draw of m trees and sigma.
    rb$yhat.test[1,] # a draw from posterior predictive distribution for each 

  # if tt == "anova"
    yimp <- apply(rb$yhat.test, 2, sample, size=1) # from each predictive distribution we 
                                                 # sample just one value
  
  # If tt==class  
    yimp_prob <- plogis(yimp)
    yimp <- ifelse(yimp_prob <= .5, 0, 1)
    
    
  Z0[wy, j] <- yimp
  
  
  
  
  
  
  
  
  
  
  
  
  
  