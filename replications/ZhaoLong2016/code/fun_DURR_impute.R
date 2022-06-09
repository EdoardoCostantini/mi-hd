### Title:    imputeHD-comp impute w/ Regu Freq Regr DURR
### Author:   Anonymized for peer review
### Created:  2020-01-8
### Modified: 2020-02-24
### Notes:    function to impute data according to DURR method following 
###           reference papers Zhao Long 2016 (for univariate miss) and 
###           Deng et al 2016 (for multivariate miss)
  
impute_DURR <- function(Xy_mis, cond, chains=5, iters=1, reg_type="lasso"){
  ## Description
  ## Packages required by function
  # library(glmnet)     # for regularized regressions
  # library(tidyverse)  # for piping
  # library(caret)      # for trainControl
  # library(PcAux)      # for iris data
  ## Other functions needed:
  # source("./init.R")
  # source("./functions.R")
  ## Input: a simulation condition, number of chians and iterations, reg type
  ## For internals:
  # cond     = conds[1,]

  # # Z      = iris2[,-c(1,7)] # alternative
  # Z0       = init_dt_i(Z, missing_type(Z)) # initialized dataset
  # reg_type = c("el", "lasso")[2]           # imputation model penality type
  # iters    = 1                             # number of iterations
  # chains   = 5                             # number of imputed datasets
  ## output: an object containing iters number of imputed datasets (imputed_datasets)
  
  ## Body
  # Prep objects
  Z        = Xy_mis
  Zm <- init_dt_i(Z, missing_type(Z))
  
  r  <- colSums(!is.na(Z)) # variable-wise count of observed values
  vmis_ind <- r != nrow(Z)
  p_imp <- sum(vmis_ind)
  p  <- ncol(Z) # number of variables [indexed with j]
  
  imputed_datasets <- vector("list", chains)
  names(imputed_datasets) <- seq(1, chains)
  
  for (m in 1:chains) {
    
    print(paste0("DURR Chain: ", m, " of ", chains))
    
    for(i in 1:iters) {
      
      # print(paste0("Iteration status: ", i, " of ", iters))
      
      for (j in 1:p_imp) { 
        # j <- 5
        J <- which(vmis_ind)[j]
        glmfam <- detect_family(Zm[, J])
        
        Zy <- process_4DURR(Z = Z, Zm = Zm, j_th = J, parms = parms)
        
        Wstar_mj <- Zy$Wstar_mj
        zstar_j <- Zy$zstar_j
        Wm_mj <- Zy$Wm_mj
        
# Fit regularized regressions ---------------------------------------------
        
        ## Lasso
        if(reg_type == "lasso"){
          regu.mod <- rr_est_lasso(X = Wstar_mj, y = zstar_j, parms = parms)
        }
        
        ## Elastic net
        if(reg_type == "el"){
          regu.mod <- rr_est_elanet(X = Wstar_mj, y = zstar_j, parms = parms)
        }
        
# Impute ------------------------------------------------------------------
        
        if(glmfam == "gaussian"){
          zm_j <- imp_gaus_DURR(regu.mod, Wstar_mj, zstar_j, Wm_mj, parms)
        }
        if(glmfam == "binomial"){
          zm_j <- imp_dich_DURR(regu.mod, Wstar_mj, zstar_j, Wm_mj, parms)
        }
        if(glmfam == "multinomial"){
          zm_j <- imp_multi_DURR(regu.mod, Wstar_mj, zstar_j, Wm_mj, parms)
        }

# Append ------------------------------------------------------------------

        Zm[is.na(Z[, J]), J] <- zm_j
        
      }

# Store m-th imputed data -------------------------------------------------

      imputed_datasets[[m]] <- Zm # Store imputed dataset for iteration
      
    }

  }
  
  return(imputed_datasets)
  
}

  