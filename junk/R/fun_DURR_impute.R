### Title:    imputeHD-comp impute w/ Regularized Frequentiest Regressions DURR
### Author:   Edoardo Costantini
### Created:  2020-JAN-8
### Modified: 2020-JAN-9
### Notes:    function to impute data according to DURR method following reference paper Deng et al 2016
  
impute_DURR <- function(Z, Z0, iters=5, reg_type="el"){
  # # Packages required by function
  #   library(glmnet)     # for regularized regressions
  #   library(tidyverse)  # for piping
  #   library(caret)      # for trainControl
  #   library(PcAux)      # for iris data
  # # Functions needed for example inputs processing (done to create the arguments of impute_DURR)
  #     source("./functions_allpurp.R")
  # # input: a dataset with missing values, an intialized version, 
  #          and some details about the procedure
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
    for (j in 1:p) { # for j-th variable w/ missing values in p number of variables w/ missing values 
                     # skipping variables without missing with the if
      # j <- 6
      if(r[j] != nrow(Z)){ # perform only for variables that have missing values
        ry <- !is.na(Z[, j]) # indexing observed values on y
        wy <- !ry            # indexing missing values on y
        Wm_j  <- Zm[,-j]
        zm_j  <- Zm[,j]
        Wm_mj <- Wm_j[wy, ]  # Xs for j-th var w/ originally missing values [OBSERVED and m-1 IMPUTATIONS]
        zm_mj <- zm_j[wy]    # imputation of j-th var at m-1 iteration
        
        ## Generate bootstrap sample
        indx_boSample <- sample(1:nrow(Wm_j), nrow(Wm_j),replace = TRUE)
        W.star.m_j    <- Wm_j[indx_boSample,]
        z.star_j      <- zm_j[indx_boSample]
        
        # but select only real observed values out of such sample
        z.star_j_obs     <- z.star_j[!is.na(Z[indx_boSample, j])]         # OBSERVED values on j-th
        W.star.m_j_obs   <- W.star.m_j[!is.na(Z[indx_boSample, j]), ] # Xs for j-th var originally observed values [OBSERVED and m-1 IMPUTATIONS]
        
        # Data for regularized models
        X <- model.matrix(z.star_j_obs~., data.frame(z.star_j_obs, W.star.m_j_obs))[,-1]
        y <- z.star_j_obs
        
        glmfam   <- detect_family(y) # regularized family type

        ## Fit regularized regression (elastic net)
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
        
        # Impute
        if(glmfam == "gaussian"){
          # x <- x[, row.names(coef(regu.mod))[as.vector(!coef(regu.mod)==0)]]
          s2.hat.m_j <- mean((predict(regu.mod, X) - z.star_j_obs)**2) # according to paper this is the estimate
          x4pred     <- model.matrix(zm_mj~., data.frame(zm_mj, Wm_mj))[,-1] # create a prediction matrix (for possible dummy coded needed)
          z.m_j_mis  <- rnorm(nrow(Wm_mj), predict(regu.mod, x4pred), sqrt(s2.hat.m_j))
        }
        if(glmfam == "binomial"){
          x4pred    <- model.matrix(zm_mj~., data.frame(zm_mj, Wm_mj))[,-1] # create a prediction matrix (for possible dummy coded needed)
          py1       <- predict(regu.mod, x4pred, type = "response")         # obtain the predicted probabilties for missing values based on their original dataset other values
          z.m_j_mis <- rbinom(nrow(py1), 1, py1)                            # sample from binomail distirbution with the drawn probabilities
          # Clumsy way of fixing the levels
          z.m_j_mis <- factor(z.m_j_mis) # return to original labels
          levels(z.m_j_mis) <- levels(y)
        }
        if(glmfam == "multinomial"){
          x4pred    <- model.matrix(zm_mj~., data.frame(zm_mj, Wm_mj))[,-1] # create a prediction matrix (for possible dummy coded needed)
          py        <- predict(regu.mod, x4pred, type = "response") # obtain the predicted probabilties for missing values based on their original dataset other values
          # rmultinom(n=1, size=1, py[1,,])
          DV_location = t(apply(py, 1, rmultinom, n = 1, size = 1)) # 1 for the category in which is most 
                                                                    # likely that an observation is
          colnames(DV_location) <- levels(y)
          z.m_j_mis <- apply(DV_location, 1, function(x) names( which(x==1) ))
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
  
  