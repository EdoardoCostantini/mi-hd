### Title:    imputeHD-comp impute w/ Regularized Frequentiest Regressions DURR
### Author:   Edoardo Costantini
### Created:  2019-NOV-27
### Modified: 2019-DEC-20
### Notes:    reference paper is Deng et al 2016

# load packages
library(tidyverse)
library(glmnet)     # for regularized regressions
library(mice)      
library(PcAux)      # for iris2 dataset
library(caret)      # for train function that you use for crossvalidation of EN parameters
library(doParallel) # to parallelize cross validation

# Prep data ---------------------------------------------------------------
# Load all purpose functions
  source("./functions_allpurp.R")
# Create using datagen function
  source("./dataGen_test.R")
    set.seed(20191120)
  dt <- missDataGen(n=100, p=10)
  dt_c <- dt[[1]] # fully observed
  dt_i <- dt[[2]] # with missings
    dim(dt_i)
    mice::md.pattern(dt_i)
  # or use iris dataset which has missing values and mix of continuous and categorical variables
  data(iris2)
  dt_i <- iris2[,2:6]
    dim(dt_i)
    mice::md.pattern(dt_i)
    row.names(dt_i)
    str(selfreport)
    selfreport$ed
    colSums(is.na(selfreport[, c("src", "age", "sex", "hr", "wr", "pop", "web", "prg", "edu")]))
    
# Match paper (see note up top) notation
  Z <- dt_i    # dataset with missing values
  p <- ncol(Z) # number of variables [INDEX with j]
  n <- nrow(Z) # number of observations
    colnames(Z) <- paste0(rep("z_", p), seq(1, p))
  r <- colSums(apply(Z, 2, function (x) {!is.na(x)} )) # vector with number of observed values in a j variable
  
  # # For a jth variable say (example notation)
  #   j = 2
  #   z_j_obs   <- Z[, j][!is.na(Z[, j])]  # observed components of j-th variable
  #   z_j_mis   <- Z[, j][is.na(Z[, j])]   # missing components of j-th variable
  #   z_mj      <- Z[, -j]                 # collection of the p − 1 variables in Z except z_j (m for minus)
  #   z_mj_obs  <- z_mj[!is.na(Z[, j]), ]  # z_mj compponents corresponding to z_j_obs
  #   z_mj_mis  <- z_mj[is.na(Z[, j]), ]   # z_mj compponents corresponding to z_j_mis
    
# Define variables with missings
  l <- missing_type(Z)$l             # number of variables needing imputation
  l_names <- missing_type(Z)$l_names # select the names of these k variables
  
# Initiliaze (by sample mean or mode category)
  set.seed(20191219)
  Z0 <- init_dt_i(Z, missing_type(Z)) # initialized dataset
  Zm <- Z0  # the dataset at Zm iteration 
            # when dataset has been itialized, m=0, so Z0 = {z0_1, z0_2, ... , z0_l, z_l+1, ... , z_p}
            # each z_j of this data will be the m-1 "previous iteration" version at the beginning of 
            # the variable loop (for j in 1:l) and the current iteration data at the end
  
# Impute dataset
  reg_type <- c("el", "las")[1]  # penality type for imputation (las not supported yet)
  iters <- 5
  imputed_datasets <- vector("list", iters)
  
  for(m in 1:iters) {
    print(paste0("Iteration status: ", m, " of ", iters))
    for (j in 1:p) { # for j-th variable w/ missing values in p number of variables w/ missing values 
                     # skipping variables without missing with the if

      if(r[j] != nrow(Z)){ # perform only for variables that have missing values
        ry <- !is.na(Z[, j]) # indexing obserrved values on y
        wy <- !ry            # indexing missing values on y
        Wm_j  <- Zm[,-j]
        zm_j  <- Zm[,j]
        Wm_mj <- Wm_j[wy, ]  # Xs for j-th var w/ originally missing values [OBSERVED and m-1 IMPUTATIONS]
        zm_mj <- zm_j[wy]    # imputation of j-th var at m-1 iteration
        
        ## Generate bootstrap sample
        indx_boSample <- sample(1:nrow(Wm_j), nrow(Wm_j),replace = TRUE)
        W.star.m_j    <- Wm_j[indx_boSample,]
        z.star_j      <- zm_j[indx_boSample]                          # not used
        
        # but select only real observed values out of such sample
        z.star_j_obs     <- zm_j[!is.na(Z[indx_boSample, j])]         # OBSERVED value son j-th
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

        ## Fit regularized regression (lasso)
        # if(reg_type == "lasso"){
        #   # Lasso Regression: choose lambda with corss validation
        #   cv.out = cv.glmnet(x[,-1], y, family = glmfam,
        #                      nfolds = 10, # as specified in paper (also default)
        #                      alpha = 1) # alpha = 1 is the lasso penality
        #   b_lambda <- cv.out$lambda.min
        # 
        #   # Fit rigde Regression with best lambda
        #   regu.mod = glmnet(x[,-1], y,
        #                     family = glmfam,
        #                     alpha  = 1,
        #                     lambda = b_lambda,
        #                     thresh = 1e-12)
        # }

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
  

# Check models ------------------------------------------------------------
  dt4check <- dt_c                  # use complete data 
  dt4check <- imputed_datasets[[5]] # use one imputed dataset
  
  # True models
  lm(y~V1+V2+V3+V4+V5, dt4check)    # linear model
  glm(y_dicho ~ V1 , data = dt4check, family = binomial) #dichotmous
  multinom(y_categ ~ V1, dt4check) #multinomial (right now I cannot really replicate this results with the lasso penality)
  
  # Lasso penalities models (change dependent variable to match interest)
  X <- model.matrix(z_4~., dt4check[,1:10])[,-1]
  y <- dt4check$z_4
  glmfam <- detect_family(y)
  
  # Lasso Regression: choose lambda with corss validation
  cv.out = cv.glmnet(x, y, family = glmfam,
                     alpha = 1) # alpha = 1 is the lasso penality
  b_lambda <- cv.out$lambda.min
  
  # Fit rigde Regression with best lambda
  lasso.mod = glmnet(x, y,
                     family = glmfam, 
                     alpha  = 1,
                     lambda = b_lambda,
                     thresh = 1e-12)
  coef(lasso.mod)
  lasso.mod$beta
  
  # Fit elastic net model
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
  coef(regu.mod)


  