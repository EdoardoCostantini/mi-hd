### Title:    imputeHD-comp impute w/ Regu Freq Regr IURR
### Author:   Anonymized for peer review
### Created:  2020-02-23
### Modified: 2020-02-23
### Notes:    function to impute data according to IURR method following 
###           reference papers Zhao Long 2016 (for univariate miss) and 
###           Deng et al 2016 (for multivariate miss). The function is a bit
###           more complex than it needs: it supports multivariate miss.
  
impute_IURR <- function(Xy_mis, cond, chains=5, reg_type="lasso"){
  # # Description
  # # Packages required by function
  # library(glmnet)     # for regularized regressions
  # library(tidyverse)  # for piping
  # library(caret)      # for trainControl
  # library(PcAux)      # for iris data
  # # Other functions needed:
  # source("./init.R")
  # source("./functions.R")
  # # Input: a simulation condition, number of chians and iterations, reg type
  # # For internals:
  # cond     = conds[9, ]
  # Xy <- genData(cond, parms)
  # # Xy_mis  = iris2[,-c(1,7)] # alternative
  # reg_type = c("el", "lasso")[2]           # imputation model penality type
  # iters    = 1                             # number of iterations
  # chains   = 5                             # number of imputed datasets
  # # output: an object containing iters number of imputed datasets (imputed_datasets)
  
  ## Body
  # Prep objects
  Z <- Xy_mis
  Zm <- Z
  p  <- ncol(Z) # number of variables [indexed with j]
  p_imp <- 1 # univariate missingnes case
  
  imputed_datasets <- vector("list", chains)
  names(imputed_datasets) <- seq(1, chains)
  
  for (m in 1:chains) {
    
    print(paste0("IURR Chain: ", m, " of ", chains))
    
    j <- which(colSums(is.na(Xy_mis)) != 0)
    O <- 1 * (!is.na(Xy_mis))
    
    z_j_obs  <- y_obs <- as.vector(Xy_mis[O[, j] == 1, j])
    zm_mj    <- y_mis <- as.vector(Xy_mis[O[, j] == 0, j])
    Wm_j_obs <- X_obs <- as.matrix(Xy_mis[O[, j] == 1, -j])
    Wm_mj    <- X_mis <- as.matrix(Xy_mis[O[, j] == 0, -j])
    X <- as.matrix(Xy_mis[, -j])
    

    # Fit regularized regression ----------------------------------------------
    # Lasso
    if(reg_type == "lasso"){
      regu.mod <- rr_est_lasso(X = Wm_j_obs, y = z_j_obs, parms = parms)
    }
    # Elastic net
    if(reg_type == "el"){
      regu.mod <- rr_est_elanet(X = Wm_j_obs, y = zstar_j, parms = parms)
    }
    
    # Impute ------------------------------------------------------------------
    zm_j <- imp_gaus_IURR(model = regu.mod,
                          X_tr = Wm_j_obs, y_tr = z_j_obs,
                          X_te = Wm_mj, y_te = zm_mj, 
                          parms = parms)

    
    # Append ------------------------------------------------------------------
    Zm[O[, j] == 0, j] <- zm_j
    
    # Store m-th imputed data -------------------------------------------------
    imputed_datasets[[m]] <- Zm # Store imputed dataset for iteration
      
  }
  
    return(imputed_datasets)
  
  }
  