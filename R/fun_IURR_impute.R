### Title:    imputeHD-comp impute w/ Regularized Frequentiest Regressions DURR
### Author:   Edoardo Costantini
### Created:  2020-JAN-8
### Modified: 2020-JAN-9
### Notes:    function to impute data according to IURR method following reference paper Deng et al 2016

# load packages

# Prep data ---------------------------------------------------------------
    
impute_IURR <- function(Z, Z0, iters=5, reg_type="el"){
  # # Packages required by function
  #   library(glmnet)    # for regularized regressions
  #   library(nnet)      # for multinomial logistic
  #   library(tidyverse) # for piping
  #   library(caret)     # for trainControl
  #   library(PcAux)     # for iris data
  # # Functions needed for example inputs processing (done to create the arguments of impute_IURR)
  #     source("./functions_allpurp.R")
  # # input: a dataset with missing values, an intialized version, 
  #        and some details about the procedure
  # #   examples:
  #       data(iris2)
  #       Z        <- iris2[,-c(1,7)]               # dataset with missing values
  #       Z0       <- init_dt_i(Z, missing_type(Z)) # initialized dataset
  #       reg_type <- c("el", "las")[1]             # penality type for imputation model (exclusive lasso not supported yet)
  #       iters    <- 5                             # number of iterations
  # # output: an object containing iters number of imputed datasets (imputed_datasets)
    
  # Body
  Zm <- Z0
  r <- colSums(apply(Z, 2, function (x) {!is.na(x)} )) # vector with number of observed values in a j variable
  p <- ncol(Z) # number of variables [indexed with j]
  imputed_datasets <- vector("list", iters)
  
  for(m in 1:iters) {
    # m <- 1
    print(paste0("Iteration status: ", m, " of ", iters))
    for (j in 1:p) {                   # for j-th variable in p number of variables
      # j <- 4
      if(r[j] != nrow(Z)){             # perform only for variables that have missing values
        ## Step 0. Prepare data for j-th variable imputation
        ry <- !is.na(Z[, j]) # indexing obserrved values on y
        wy <- !ry            # indexing missing values on y
        Wm_j  <- Zm[,-j]     # predictors for imp model from inizialized dataset (m-1) [ALL CASES]  
        zm_j  <- Zm[,j]      # outcome for imp model from inizialized dataset (m-1)    [ALL CASES]
        Wm_mj <- Wm_j[wy, ]  # predictor rows for cases with missing z_j               [OBSERVED and IMPUTED]
        zm_mj <- zm_j[wy]    # current (m) draw of z_j for cases with missing z_j      [IMPUTED/MISSING CASES]
        
        # Select cases for imputation model
        z_j_obs  <- Zm[, j][!is.na(Z[, j])]   # observed components of j-th variable [OBSERVED VALUES]
        Wm_j_obs <- Zm[, -j][!is.na(Z[, j]),] # current (m) components of the predictors cooresponding 
                                              # to observed cases on j-th variable [OBSERVED and IMPUTED VALUES]
        
        ## Active set selection thorugh regularized regression (Elastic Net penality) ##
        X <- model.matrix(z_j_obs~., data.frame(z_j_obs, Wm_j_obs))[,-1]
        y <- z_j_obs
        glmfam   <- detect_family(y) # regularized family type
        
        # Fit regularized regression (elastic net)
        if(reg_type == "el"){
          # Set training control
          train_control <- trainControl(method = "cv",
                                        number = 10,
                                        # repeats = 5,
                                        search = "random", # NEEDS CHECKING: what is this search method?
                                        verboseIter = T)
          
          # CV
          el_cv <- train(y ~., 
                         data = cbind(y, X), 
                         method = "glmnet",
                         family = glmfam, 
                         type.multinomial = "grouped",
                         trControl = train_control,
                         preProcess = c("center", "scale"),
                         tuneLength = 25
          )
          
          # Train model
          regu.mod <- glmnet(X, y,
                             family = glmfam,
                             type.multinomial = "grouped", # if glmfam is multinomial, otherwise does not apply
                             alpha = el_cv$bestTune[1],
                             lambda = el_cv$bestTune[2])
        }

        # Conituous and dichotomous case
        if(glmfam != "multinomial"){
          # Define Active Set for binomial and gaussian distributed variables under imputation
          lasso.coef <- as.matrix(coef(regu.mod))
          active.set <- row.names(lasso.coef)[lasso.coef != 0][-1]
          W_Sjm_obs  <- cbind(1,X[, active.set]) # reusing x object without making another model matrix
                                                 # the unclusion of the first column of 1s is for the 
                                                 # following mle fit done manually
          W_Sjm_mis  <- cbind(1,model.matrix(zm_mj~., data.frame(zm_mj, Wm_mj))[, active.set])
          
          if(glmfam == "gaussian"){
            # Fit ML regression
            # Obtain coefficeints estiamtes (to faciliate the second step)
            lm_fit <- lm(y ~ W_Sjm_obs[,-1]) # lm does not want the desing matrix with 1 for intercepts
            startV <- c(coef(lm_fit), sigma(lm_fit))
            # Obtain Hessian matrix (for standard error of the error variance)
            ml.fit <- optim(startV, # initial values (use the guesses from a lm fit to make it always ok)
                                    # what matters is simply obtaining the standard error of the error
                                    # variance that lm does not provide. Hence, I can make it easy for the pc
                                    # and provide the parameters value as close to reality as possible and 
                                    # just exploit what I need of this procedure.
                            unigaus.lf.dnorm,          # costum log likelihood function (version w/ dnorm)
                            method="BFGS",
                            hessian=T,
                            y = y, X = W_Sjm_obs)

            theta.hat.m_MLE <- ml.fit$par
              # vector of estiamted reg coef at this iteration
              # slight variation from paper notation is that I'm including the error varaince parameter in this
              # vector theta although in the paper theta are just slopes and intercepts
            
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
            # Clumsy way of fixing the levels
            z.m_j_mis <- factor(z.m_j_mis) # return to original labels
            levels(z.m_j_mis) <- levels(y)
          }
          
        }
        
        # Multinomial case
        if(glmfam == "multinomial"){
          # Define Active Set
          lasso.coef <- as.matrix(coef(regu.mod)[[1]])
          active.set <- row.names(lasso.coef)[lasso.coef != 0][-1] # exclude intercept missing name
          W_Sjm_obs  <- model.matrix(z_j_obs~., data.frame(z_j_obs, Wm_j_obs))[, c("(Intercept)", active.set)]
          W_Sjm_mis  <- model.matrix(zm_mj~., data.frame(zm_mj, Wm_mj))[, c("(Intercept)", active.set)]

          MLfit <- multinom(y ~ 0 + W_Sjm_obs) # first Y category as reference
          
          theta.hat.m_MLE <- as.vector(t(coef(MLfit)))
          Sigma.hat.m.MLE <- vcov(MLfit)
          
          # Sample parameters for prediction/imputation
          theta.hat.m_j <- MASS::mvrnorm(1, theta.hat.m_MLE, Sigma.hat.m.MLE)
          
          coef_matrix <- matrix(theta.hat.m_j, ncol = ncol(coef(MLfit)), byrow = T)
            # each column is for a different parameter, each row for a different category
          
          # Get imputations
          K <- length(unique(y)) # number of categories of variable under imputation
          logit <- matrix(NA, nrow = length(zm_mj), ncol = K-1)
          for (k in 1:(K-1)) {
            if(ncol(coef_matrix)==1){
              logit[,k] <- W_Sjm_mis * coef_matrix[k,]
            } else {
              logit[,k] <- W_Sjm_mis %*% coef_matrix[k,]
            }
          }
          
          # Augmente matrix of logits to obtain probabilities for all categories directly
          logit_A <- cbind(rep(0, nrow(as.data.frame(W_Sjm_mis))), logit)
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
  return(imputed_datasets)
}

  