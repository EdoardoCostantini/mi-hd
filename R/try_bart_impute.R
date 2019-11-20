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
    
  # Create using datagen function
  source("./dataGen_test.R")
  set.seed(20191120)
  dt <- missDataGen(n=100, p=500)
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

# Sequential BART
  iters <- 10
  fast_diagnostic <- matrix(rep(NA, iters*7), ncol = 7)
for(i in 1:iters) {
  print(paste0("Iteration status: ", i, " of ", iters))
  for (k in 1:K) {
    print(paste0("Variable under imputation: ", K_names[k], " (", k, " of ", K, ")" ))
    # Prep Data
    indx_k <- names(dt_aug) %in% K_names[k]               # to index x_k
    ry <- !is.na(dt_i[, K_names[k]])                      # response variable T = observed, F = NA
    data4tree <- arrange_dt_x_k(dt_aug, K_names[k])[ry, ] # arrange and select Xobs (values of Xno_k for which X_k is observed)
    Xobs <- data4tree[, -1]                               # after rearranging, the variable under imputation is always first
      dim(Xobs)
    yobs <- data4tree[, 1]
      length(yobs)
    Xmiss <- arrange_dt_x_k(dt_aug, K_names[k])[!ry, -1]     # observed and previously imputed / initialized values on all other variables
      dim(Xmiss)
    # Define tree type (classification or regression) based on nature of variable under imputation
    #tt <- tree_type(data4tree$x_k)
    
    # Grow tree
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
    dt_aug[row.names(Xmiss), indx_k] <- bt$yhat.test.mean
  }
   new_impute <- dt_aug %>% filter(is.na(dt_i$yinc)) %>% select(1:7)
   original <- dt_c %>% filter(is.na(dt_i$yinc)) %>% select(1:7)
   fast_diagnostic[i,] <- colSums(new_impute - original)
}
  
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
  
  
  