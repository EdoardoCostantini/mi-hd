### Title:    imputeHD-comp impute w/ Regularized Frequentiest Regressions DURR
### Author:   Edoardo Costantini
### Created:  2019-DEC-2
### Modified: 2019-DEC-4
### Notes:    reference paper is Deng et al 2016

# load packages
library(tidyverse)
library(dplyr)
library(glmnet)    # for regularized regressions
library(nnet)      # for multinomial logistic
library(mice)      
library(caret)
library(doParallel)

# Prep data ---------------------------------------------------------------
# Load all purpose functions
  source("./functions_allpurp.R")
# Create using datagen function
  source("./dataGen_test.R")
    set.seed(20191120)
  dt <- missDataGen(n=500, p=500)
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
  print("Data initialization")
  
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
  iters <- 5
  imputed_datasets <- vector("list", iters)
  for(m in 1:iters) {
    print(paste0("Iteration status: ", m, " of ", iters))
    for (j in 1:p) {                   # for j-th variable in p number of variables
      if(r[j] != nrow(Z)){             # perform only for variables that have missing values
        j <- 4
        ## Step 0. Prepare data for j-th variable imputation
        Wm_j  <- Zm[,-j]               # predictors for imp model from inizialized dataset (m-1) [ALL CASES]
        zm_j  <- Zm[,j]                # outcome for imp model from inizialized dataset (m-1)    [ALL CASES]
        Wm_mj <- Wm_j[is.na(Z[, j]), ] # predictor rows for cases with missing z_j               [MISSING CASES]
        zm_mj <- zm_j[is.na(Z[, j])]   # current (m) draw of z_j for cases with missing z_j      [MISSING CASES]
        
        # Select cases for imputation model
        z_j_obs  <- Zm[, j][!is.na(Z[, j])]   # observed components of j-th variable [OBSERVED VALUES]
        Wm_j_obs <- Zm[, -j][!is.na(Z[, j]),] # current (m) components of the predictors cooresponding 
                                              # to observed cases on j-th variable [OBSERVED and IMPUTED VALUES]
        
        ## Active set selection thorugh regularized regression (Elastic Net penality) ##
        x <- model.matrix(z_j_obs~., data.frame(z_j_obs, Wm_j_obs))
        y <- z_j_obs
        glmfam <- detect_family(y)
        
        # Cross validation for alpha and lambda
        cl <- makePSOCKcluster(5) # Set up parallel computing (#cores)
        registerDoParallel(cl)
        
        # Tuning parameters range
        # alpha.seq <- seq(0,1, .1)    # choosing .1 and .9 as bound to force difference from ridge and lasso
        # lambda.seq <- seq(0,10, .05)   # choosing 10 as upper bound because
        # desired.grid <- expand.grid(alpha.seq, lambda.seq)
        #   colnames(desired.grid) <- c("alpha", "lambda")
          
        # Cross-validate
        model <- train(
          y ~., data = data.frame(y, x[,-1]), method = "glmnet",
          family = glmfam, type.multinomial = "grouped", # type.multinomial is used only if glmfam = "multinomial"
          trControl = trainControl("cv", number = 10),   # 10-fold corssvalidation
          #tuneGrid = desired.grid
          tuneLength = 25 # how many alpha values tried between .1 and 1, and how many lambda values
        )
        stopCluster(cl) # terminate parallel computing
        best_al <- model$bestTune
        
        # Fit elastic net with corssvalidated values (repetation)
        el.mod <- glmnet(x[,-1], y,
                        family = glmfam,
                        type.multinomial = "grouped",
                        alpha  = best_al["alpha"],
                        lambda = best_al["lambda"],
                        thresh = 1e-12)

        # Define Active Set for binomial and gaussian distributed variables under imputation
        if(glmfam != "multinomial"){ # multinomial not currently supported
          lasso.coef <- as.matrix(coef(el.mod))
          active.set <- row.names(lasso.coef)[lasso.coef != 0]
          W_Sjm_obs  <- x[, active.set] # reusing x object without making another model matrix
          W_Sjm_mis  <- model.matrix(zm_mj~., data.frame(zm_mj, Wm_mj))[, active.set]
        }
        if(glmfam == "gaussian"){
          ## Fit ML regressions
          # Manual Way
          ml.fit <- optim(rep(1, ncol(W_Sjm_obs)+1), # initial values ("1" for each, +1 for the error variance)
                     unigaus.lf,                     # costum log likelihood function
                     #method="BFGS", 
                     hessian=T,
                     y = y, X = W_Sjm_obs)
          
          theta.hat.m_MLE <- ml.fit$par 
            # vector of estiamted reg coef at this iteration
            # slight variation from paper notation is that I'm including the error varaince parameter in this
            # vector theta although in the paper theta are just sloeps and intercepts
          Sigma.hat.m.MLE <- solve(ml.fit$hessian) # variance covariance matrix of the estiamted parameters
          
          # Sample parameters for prediction/imputation
          theta.hat.m_j <- MASS::mvrnorm(1, theta.hat.m_MLE, Sigma.hat.m.MLE)
          # Get imputations
          z.m_j_mis <- rnorm(n = nrow(W_Sjm_mis),
                            mean = W_Sjm_mis %*% theta.hat.m_j[1:(length(theta.hat.m_j)-1)],
                            sd = sqrt(tail(theta.hat.m_j,1)) )
        }
        
        if(glmfam == "binomial"){
          MLfit <- glm(y ~ 0 + W_Sjm_obs, family = glmfam) # 0+ because X is already desing matrix (p + intercept)
          theta.hat.m_MLE <- as.vector(MLfit$coefficients)
          Sigma.hat.m.MLE <- vcov(MLfit)
          
          # Sample parameters for prediction/imputation
          theta.hat.m_j <- MASS::mvrnorm(1, theta.hat.m_MLE, Sigma.hat.m.MLE)
          
          # Get imputations
          lin_term  <- W_Sjm_mis %*% theta.hat.m_j # obtain the predicted probabilties for missing values based on their original dataset other values
          z.m_j_mis <- rbinom(n = length(lin_term), size=1, 
                              prob=exp(lin_term)/(1+exp(lin_term)))
          z.m_j_mis <- factor(z.m_j_mis, labels = levels(y)) # return to original labels
        }
        
        if(glmfam == "multinomial"){
          # Define Active Set
          lasso.coef <- as.matrix(coef(el.mod)[[1]])
          active.set <- row.names(lasso.coef)[lasso.coef != 0][-1] # exclude intercept missing name
          W_Sjm_obs  <- x[, c("(Intercept)", active.set)] # reusing x object without making another model matrix
          W_Sjm_mis  <- model.matrix(zm_mj~., data.frame(zm_mj, Wm_mj))[, c("(Intercept)", active.set)]
          
          MLfit <- multinom(y ~ 0 + W_Sjm_obs) # first Y category as reference
          
          theta.hat.m_MLE <- as.vector(t(coef(MLfit)))
          Sigma.hat.m.MLE <- vcov(MLfit)
          
          # Sample parameters for prediction/imputation
          theta.hat.m_j <- MASS::mvrnorm(1, theta.hat.m_MLE, Sigma.hat.m.MLE)
          
          coef_matrix <- matrix(theta.hat.m_j, ncol = ncol(coef(MLfit)), byrow = T)
          
          # Get imputations
          K <- length(unique(y)) # number of categories of variable under imputation
          logit <- matrix(NA, nrow = length(zm_mj), ncol = K-1)
          for (k in 1:(K-1)) {
            logit[,k] <- W_Sjm_mis %*% coef_matrix[k,]
          }
          # Augmente matrix of logits to obtain probabilities for all categories directly
          logit_A <- cbind(rep(0, nrow(W_Sjm_mis)), logit)
          oddsratio_A <- exp(logit_A)
          vProb <- oddsratio_A/(rowSums(oddsratio_A))
          # Predict missing values
          DV_location = t(apply(vProb, 1, rmultinom, n = 1, size = 1))
          DV <- apply(DV_location, 1, function(x) which(x==1))
          z.m_j_mis <- factor(DV, labels = levels(y))
        }
        # Append
        Zm[is.na(Z[, j]), j] <- z.m_j_mis
      }
    }
    
    # Store imputed dataset at this iteration
    imputed_datasets[[m]] <- Zm
  }

  
# Model checks on complete data
# Original data
  dt_c
  y <- dt_c$y
  x <- model.matrix(y~., dt_c)
  head(x)
# True linear reg model
  lm(y~x[,c("V1","V2","V3","V4","V5")])

# Fit elastic net
  # Cross validation for alpha and lambda
  cl <- makePSOCKcluster(5) # Set up parallel computing (#cores)
  registerDoParallel(cl)

  # Cross-validate
  model <- train(
    y ~., data = data.frame(y, x[,-1]), method = "glmnet",
    family = glmfam, type.multinomial = "grouped", # type.multinomial is used only if glmfam = "multinomial"
    trControl = trainControl("cv", number = 10),   # 10-fold corssvalidation
    #tuneGrid = desired.grid
    tuneLength = 25 # how many alpha values tried between .1 and 1, and how many lambda values
  )
  stopCluster(cl) # terminate parallel computing
  best_al <- model$bestTune
  
  # Fit elastic net with corssvalidated values (repetation)
  el.mod <- glmnet(x[,-1], y,
                   family = glmfam,
                   type.multinomial = "grouped",
                   alpha  = best_al["alpha"],
                   lambda = best_al["lambda"],
                   thresh = 1e-12)
  coef(el.mod)
  md.pattern(model$results)

  model$results[is.na(model$results$RsquaredSD), ]
  