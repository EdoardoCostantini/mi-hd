### Title:    imputeHD-comp impute w/ Regularized Frequentiest Regressions DURR
### Author:   Edoardo Costantini
### Created:  2019-NOV-27
### Modified: 2019-DEC-2
### Notes:    reference paper is Deng et al 2016

# load packages
library(tidyverse)
library(dplyr)
library(mice)      
library(monomvn)     # For bayesian lasso
library(EBEN)        # authored by one of the co-authors to the papers on logistic Blasso
library(EBglmnet)    # Logistic bayesian lasso (functions do not provide posterior distributions)

# Prep data ---------------------------------------------------------------
# Load all purpose functions
  source("./functions_allpurp.R")
# Create using datagen function
  source("./dataGen_test.R")
    set.seed(20191120)
  dt <- missDataGen(n=1e3, p=10)
  dt_c <- dt[[1]] # fully observed
  dt_i <- dt[[2]] # with missings
    dim(dt_i)
    mice::md.pattern(dt_i)
    
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
  
# Imputed dataset
  iters <- 5 # iterate until convergence (not clear what convergence is)
  imputed_datasets <- vector("list", iters)
  for(m in 1:iters) {
    print(paste0("Iteration status: ", m, " of ", iters))
    for (j in 1:p) { # for j-th variable w/ missing values in p number of variables w/ missing values 
                     # skipping variables without missing with the if
      if(r[j] != nrow(Z)){
        # perform only for variables that have missing values
        # iv and dv are defined in terms of the imptuation model: dv is the variable under imputation
        Wm_j  <- Zm[,-j] # predictors full data (observed + imputed)
        zm_j  <- Zm[,j] # dv full data (observed + imputed)
        
        Wo_mj <- Wm_j[!is.na(Z[, j]), ] # predictor values for those observed on z_j
        Wm_mj <- Wm_j[is.na(Z[, j]), ]  # mj = -j, predictors values for those missing on z_j
        
        zo_mj <- zm_j[!is.na(Z[, j])]   # predictor values for those observed on z_j
        zm_mj <- zm_j[is.na(Z[, j])]   # current (m) draw of z_j for cases with missing z_j [MISSING CASES]
        
        ### NEW ###
        
        Wo_mj_x <- model.matrix(zo_mj~., data.frame(zo_mj, Wo_mj))[,-1]
        y <- zo_mj
        glmfam <- detect_family(y)
        
        if(glmfam == "gaussian"){
          lasso.mod <- blasso(X = Wo_mj_x, y = y,
                              T=1000,      # I'm not sure that the function does burn-in itself
                              thin=NULL,   # automatic skip before sample is collected,
                              lambda2 = 0, # needs tp be estiamted with maximum likelihood principle CHECK
                              # sigma**2 prior: IG(a,b)
                              ab = c(.1,.1), # alpha (shape) parameter and the beta (scale) parameter for the 
                                           # IG(a,b) for the variance parameter s2
                              # lamba prior: Gamma(r,s)
                              rd = c(.01, .01), # the alpha (shape) parameter and beta (rate) parameter 
                                                # G(r,delta) for the lambda2 parameter under the lasso model
                              # rho prior: Beta(g,h)
                              mprior = c(1,1),  # Bin(m|n=M,p) prior for number of 0 coefficients, where p~Beta(g,h)
      
                              verb=1
          )

          # CHECK parameters match what advised in ZhaoLong2013
          sample_index <- sample(1:nrow(lasso.mod$beta), 1)
          theta.hat.m_j <- lasso.mod$beta[sample_index, ] # random draw from posterior distribution of theta.hat.m_j
          
          Wm_mj_x <- model.matrix(zm_mj~., data.frame(zm_mj, Wm_mj))[,-1]
          z.m_j_mis <- Wm_mj_x %*% theta.hat.m_j + rnorm(1, 0, sqrt(lasso.mod$s2[sample_index])) # sample from predictive distribution
        }
        if(glmfam == "binomial"){
          err1 <- paste0(colnames(dt_i)[j], " is a dichotomous variable. Currently not supported: returning initial guess")
          print(err1)
          z.m_j_mis <- Zm[is.na(Z[, j]), j]
        }
        if(glmfam == "multinomial"){
          err1 <- paste0(colnames(dt_i)[j], " is a polytomous variable. Currently not supported: returning initial guess")
          print(err1)
          z.m_j_mis <- Zm[is.na(Z[, j]), j]
        }
        
        # Append
        Zm[is.na(Z[, j]), j] <- z.m_j_mis
      }
    }
      
    # Store imputed dataset at this iteration
    imputed_datasets[[m]] <- Zm
  }
  