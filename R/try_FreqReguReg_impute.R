
### Title:    imputeHD-comp impute w/ Regularized Frequentiest Regressions
### Author:   Edoardo Costantini
### Created:  2019-NOV-27
### Modified: 2019-NOV-27

# load packages
library(tidyverse)
library(dplyr)
library(glmnet)    # for regularized regressions

# Prep data ---------------------------------------------------------------

# Create using datagen function
  source("./dataGen_test.R")
  set.seed(20191120)
  dt <- missDataGen(n=20, p=6)
  dt_c <- dt[[1]] # fully observed
  dt_i <- dt[[2]] # with missings
    dim(dt_i)
    mice::md.pattern(dt_i)
  # or use iris dataset which has missing values and mix of continuous and categorical variables
  dt_i <- iris2
    dim(dt_i)
    mice::md.pattern(dt_i)
    row.names(dt_i)
    
# Match paper notation
  Z <- dt_i    # dataset with missing values
  p <- ncol(Z) # number of variables [INDEX with j]
  n <- nrow(Z) # number of observations
    colnames(Z) <- paste0(rep("z_", p), seq(1, p))
  r <- colSums(apply(Z, 2, function (x) {!is.na(x)} )) # vector with number of observed values in a j variable
  
  # Define variables with missings
  l <- ncol(Z)-sum(tail(mice::md.pattern(Z),1) == 0) # number of variables needing imputation
  l_names <- names(which(colSums(apply(Z, 2, is.na)) != 0)) # select the names of these k variables
  iters <- 5
  
  # For a jth variable say (example)
    j = 2
    z_j_obs   <- Z[, j][!is.na(Z[, j])]  # observed components of j-th variable
    z_j_mis   <- Z[, j][is.na(Z[, j])]   # missing components of j-th variable
    z_mj      <- Z[, -j]                 # collection of the p − 1 variables in Z except z_j (m for minus)
    z_mj_obs  <- z_mj[!is.na(Z[, j]), ]  # z_mj compponents corresponding to z_j_obs
    z_mj_mis  <- z_mj[is.na(Z[, j]), ]   # z_mj compponents corresponding to z_j_mis
            
                  with(mtcars, mpg[cyl == 8  &  disp > 350])
                  # is the same as, but nicer than
                  mtcars$mpg[mtcars$cyl == 8  &  mtcars$disp > 350]
  
# Initiliaze by sample mean of variable with missing values
                  
  print("Data initialization by sample mean")
  sample_means <- apply(Z[, names(Z) %in% l_names], 2, mean, na.rm = TRUE)
  Z0 <- Z
  for (j in 1:l) {
    Z0 <- Z0 %>% dplyr::mutate_at(vars(l_names[j]), ~replace(., is.na(.), sample_means[j]))
  }
  imp <- mice::mice(Z, method = "pmm", m = 1)
  Z0 <- round(complete(imp), 2)
  Zm <- Z0  # the dataset at Zm iteration 
            # when dataset hass been itialized, m=0, so Z0 = {z0_1, z0_2, ... , z0_l, z_l+1, ... , z_p}
            # each z_j of this data will be the m-1 "previous iteration" version at the beginning of 
            # the variable loop (for j in 1:l) and the current iteration data at the end
  imputed_datasets <- vector("list", iters)
  for(m in 1:iters) {
    print(paste0("Iteration status: ", m, " of ", iters))
    for (j in 1:p) { # for j-th variable w/ missing values in p number of variables w/ missing values 
                     # skipping variables without missing with the if
      if(r[j] != nrow(Z)){ # perform only for variables that have missing values
        Wm_j  <- Zm[,-j]
        z_j   <- as.data.frame(Zm[,j])
        Wm_mj <- as.matrix(Wm_j[is.na(Z[, j]), ])
        
        ## Generate bootstrap sample
        indx_boSample <- sample(1:nrow(Wm_j), nrow(Wm_j),replace = TRUE)
        W.star.m_j   <- as.matrix(Wm_j[indx_boSample, ])
        z.star_j     <- z_j[indx_boSample, ]
        
        # but select only real observed values out of such sample
        z.star_j_obs     <- z_j[!is.na(Z[indx_boSample, j]), ]
        W.star.m_j_obs   <- W.star.m_j[!is.na(Z[indx_boSample, j]), ]
        
        ## Fit regularized regression
        # DV_type
        # RGU_type
  
        # Ridge Regression: choose lambda with corss validation
        cv.out = cv.glmnet(W.star.m_j_obs, z.star_j_obs, alpha = 0)
        b_lambda <- cv.out$lambda.min
        
        # Fit rigde Regression with best lambda
        ridge.mod = glmnet(W.star.m_j_obs, z.star_j_obs,
                           alpha  = 0, 
                           lambda = b_lambda,
                           thresh = 1e-12)
        coef(ridge.mod)
        s2.hat.m_j <- mean((predict(ridge.mod, W.star.m_j_obs) - z.star_j_obs)**2) # according to paper this is the estimate
        
        # Impute
        z.m_j_mis <- rnorm(nrow(Wm_mj), predict(ridge.mod, Wm_mj), sqrt(s2.hat.m_j))
        
        # Append
        Zm[is.na(Z[, j]), j] <- z.m_j_mis
      }
    }
    imputed_datasets[[m]] <- Zm
    #print(dt_aug %>% filter(is.na(dt_i$yinc)) %>% select(1:7))
  }
  tail(imputed_datasets, 2)