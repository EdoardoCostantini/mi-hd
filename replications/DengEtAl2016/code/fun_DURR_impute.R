### Title:    imputeHD-comp impute w/ Regu Freq Regr DURR
### Author:   Edoardo Costantini
### Created:  2020-01-8
### Modified: 2020-02-24
### Notes:    function to impute data according to DURR method following 
###           reference papers Zhao Long 2016 (for univariate miss) and 
###           Deng et al 2016 (for multivariate miss)

impute_DURR <- function(Xy_mis, chains=5, iters=5, reg_type="lasso"){
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
  # iters    = 5                             # number of iterations
  # chains        = 2                             # number of imputed datasets (how many to keep)
  ## output: an object containing iters number of imputed datasets (imputed_datasets)
  
  ## Body
  Z  <- Xy_mis
  p  <- ncol(Z) # number of variables [indexed with j]
  
  r  <- colSums(!is.na(Z)) # variable-wise count of observed values
  vmis_ind <- r != nrow(Z)
  p_imp <- sum(vmis_ind)
  
  # To store imputed values and check convergence
  imputed_values <- vector("list", p_imp)
    names(imputed_values) <- names(which(r != nrow(Z)))
  for (v in 1:p_imp) {
    imputed_values[[v]] <- matrix(rep(NA, (iters*(nrow(Z)-r[v]))), 
                                  ncol = iters)
  }
  
  # To store multiply imputed datasets
  imputed_datasets <- vector("list", iters)
    names(imputed_datasets) <- seq(1, iters)
    
  # Time performance
  start.time <- Sys.time()
  
  imp_res <- lapply(1:chains, function(x){
    Zm <- init_dt_i(Z, missing_type(Z)) # reinitialize data
    for(i in 1:iters) {
      print(paste0("DURR - Chain: ", x, "/", chains, "; Iter: ", i, "/", iters))
      for (j in 1:p_imp) {
        J <- which(vmis_ind)[j]
        glmfam <- detect_family(Zm[, J])
        
        # 1. Generate bootstrap data set (for each j!)
        Zy <- process_4DURR(Z = Z, Zm = Zm, j_th = J, parms = parms)
        
        Wstar_mj <- Zy$Wstar_mj
        zstar_j <- Zy$zstar_j
        Wm_mj <- Zy$Wm_mj
        
        # 2. Fit regularized regression (i.e. obtain estiamtes)
        ## Lasso
        if(reg_type == "lasso"){
          regu.mod <- rr_est_lasso(X = Wstar_mj, y = zstar_j, parms = parms)
        }
        
        ## Elastic net
        if(reg_type == "el"){
          regu.mod <- rr_est_elanet(X = Wstar_mj, y = zstar_j, parms = parms)
        }
        
        # Predict z_j_mis (i.e. obtain imputations)
        if(glmfam == "gaussian"){
          zm_j <- imp_gaus_DURR(regu.mod, Wstar_mj, zstar_j, Wm_mj, parms)
        }
        if(glmfam == "binomial"){
          zm_j <- imp_dich_DURR(regu.mod, Wstar_mj, zstar_j, Wm_mj, parms)
        }
        if(glmfam == "multinomial"){
          zm_j <- imp_multi_DURR(regu.mod, Wstar_mj, zstar_j, Wm_mj, parms)
        }
        
        # Append imputation
        imputed_values[[j]][, i] <- zm_j
        Zm[is.na(Z[, J]), J] <- zm_j
      }
    }  
    
    
    # Store results
    return(list(imp_dat = Zm,
                imp_val = imputed_values)
           )
  })
  
  end.time <- Sys.time()
  
  return(list(imp_res = imp_res,
              time = difftime(end.time, 
                              start.time, 
                              units = "mins"))
         )
  
}
  
  