### Title:    imputeHD-comp impute w/ Regularized Frequentiest Regressions DURR
### Author:   Edoardo Costantini
### Created:  2019-NOV-27
### Modified: 2019-DEC-2
### Notes:    reference paper is Deng et al 2016

# load packages
library(tidyverse)
library(dplyr)
library(glmnet)    # for regularized regressions
library(mice)      
library(PcAux)     # for iris2 dataset

# Prep data ---------------------------------------------------------------
# Load all purpose functions
  source("./functions_allpurp.R")
# Create using datagen function
  source("./dataGen_test.R")
    set.seed(20191120)
  dt <- missDataGen(n=300, p=100)
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
  iters <- 5
  imputed_datasets <- vector("list", iters)
  for(m in 1:iters) {
    print(paste0("Iteration status: ", m, " of ", iters))
    for (j in 1:p) { # for j-th variable w/ missing values in p number of variables w/ missing values 
                     # skipping variables without missing with the if
      if(r[j] != nrow(Z)){ # perform only for variables that have missing values
        Wm_j  <- Zm[,-j]
        zm_j   <- Zm[,j]
        Wm_mj <- Wm_j[is.na(Z[, j]), ]
        zm_mj <- zm_j[is.na(Z[, j])]
        
        ## Generate bootstrap sample
        indx_boSample <- sample(1:nrow(Wm_j), nrow(Wm_j),replace = TRUE)
        W.star.m_j    <- Wm_j[indx_boSample, ]
        z.star_j      <- zm_j[indx_boSample]
        
        # but select only real observed values out of such sample
        z.star_j_obs     <- zm_j[!is.na(Z[indx_boSample, j])]
        W.star.m_j_obs   <- W.star.m_j[!is.na(Z[indx_boSample, j]), ]
        
        ## Fit regularized regression
        x <- model.matrix(z.star_j_obs~., data.frame(z.star_j_obs, W.star.m_j_obs))[,-1]
        y <- z.star_j_obs
        glmfam <- detect_family(y)
        
        # Lasso Regression: choose lambda with corss validation
        cv.out = cv.glmnet(x, y, family = glmfam,
                           nfolds = 10, # as specified in paper (also default)
                           alpha = 1) # alpha = 1 is the lasso penality
        b_lambda <- cv.out$lambda.min
        
        # Fit rigde Regression with best lambda
        lasso.mod = glmnet(x, y,
                           family = glmfam, 
                           alpha  = 1,
                           lambda = b_lambda,
                           thresh = 1e-12)
        # lasso.coef.all <- as.data.frame(as.matrix(coef(lasso.mod)))
        # lasso.coef.sel <- data.frame(varn = row.names(lasso.coef)[lasso.coef$s0 != 0],
        #                              coef = lasso.coef[lasso.coef$s0 != 0, ])
        
        # Impute
        if(glmfam == "gaussian"){
          s2.hat.m_j <- mean((predict(lasso.mod, x) - z.star_j_obs)**2) # according to paper this is the estimate
          x4pred     <- model.matrix(zm_mj~., data.frame(zm_mj, Wm_mj))[,-1] # create a prediction matrix (for possible dummy coded needed)
          z.m_j_mis  <- rnorm(nrow(Wm_mj), predict(lasso.mod, x4pred), sqrt(s2.hat.m_j))
        }
        if(glmfam == "binomial"){
          x4pred    <- model.matrix(zm_mj~., data.frame(zm_mj, Wm_mj))[,-1] # create a prediction matrix (for possible dummy coded needed)
          py1       <- predict(lasso.mod, x4pred, type = "response") # obtain the predicted probabilties for missing values based on their original dataset other values
          z.m_j_mis <- rbinom(nrow(py1), 1, py1)  # sample from binomail distirbution with the drawn probabilities
          z.m_j_mis <- factor(z.m_j_mis, labels = levels(y)) # return to original labels
        }
        if(glmfam == "multinomial"){
          x4pred    <- model.matrix(zm_mj~., data.frame(zm_mj, Wm_mj))[,-1] # create a prediction matrix (for possible dummy coded needed)
          py        <- predict(lasso.mod, x4pred, type = "response") # obtain the predicted probabilties for missing values based on their original dataset other values
          rmultinom(n=1, size=1, py[1,,])
          DV_location = t(apply(py, 1, rmultinom, n = 1, size = 1)) # 1 for the category in which is most 
                                                                    # likely that an observation is
          colnames(DV_location) <- levels(y)
          z.m_j_mis <- apply(DV_location, 1, function(x) names( which(x==1) ))
        }
        # Append
        Zm[is.na(Z[, j]), j] <- z.m_j_mis
      }
    }
    # Print the coefficient for the linear model to simple check interative changes (short life)
      x <- model.matrix(z_4~., Zm)[,-1]
      y <- Zm$z_4
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
      lasso.coef <- as.data.frame(as.matrix(coef(lasso.mod)))
      lasso.coef.sel <- data.frame(varn = row.names(lasso.coef)[lasso.coef$s0 != 0],
                                   coef = lasso.coef[lasso.coef$s0 != 0, ])
      print(lasso.coef.sel)
      
    # Store imputed dataset at this iteration
    imputed_datasets[[m]] <- Zm
  }

# Check models ------------------------------------------------------------
  dt4check <- dt_c
  # Now using complete data, but can easly change it
  
  # True models
  lm(y~V1+V2+V3+V4+V5, dt4check) #linear model
  glm(y_dicho ~ V1 , data = dt4check, family = binomial) #dichotmous
  multinom(y_categ ~ V1, dt4check) #multinomial (right now I cannot really replicate this results with the lasso penality)
  
  # Lasso penalities models (change dependent variable to match interest)
  x <- model.matrix(y_categ~., dt4check[,1:10])[,-1]
  y <- dt4check$y_categ
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


  