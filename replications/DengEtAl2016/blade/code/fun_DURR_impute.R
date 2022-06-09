### Title:    imputeHD-comp impute w/ Regu Freq Regr DURR
### Author:   Anonymized for peer review
### Created:  2020-01-8
### Modified: 2020-02-24
### Notes:    function to impute data according to DURR method following 
###           reference papers Zhao Long 2016 (for univariate miss) and 
###           Deng et al 2016 (for multivariate miss)

impute_DURR <- function(Z, O, cond, reg_type="lasso", parms){
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
  # chains   = 2                             # number of imputed datasets (how many to keep)
  ## output: an object containing iters number of imputed datasets (imputed_datasets)
  
  # Z = Xy_mis
  # O <- as.data.frame(!is.na(Xy_mis))            # matrix index of observed values
  # reg_type = c("el", "lasso")[2] # imputation model penality type
  
  ## Body
  p  <- ncol(Z) # number of variables [indexed with j]
  
  p_imp <- length(parms$z_m_id) # variables with missing values
  r     <- colSums(O[parms$z_m_id]) # variable-wise count of observed values
  nr    <- colSums(!O[parms$z_m_id]) # variable-wise count of observed values
  
  # To store imputed values and check convergence
  imp_DURR_val <- vector("list", parms$chains)
   names(imp_DURR_val) <- seq(1:parms$chains)
    
  # Time performance
  start.time <- Sys.time()
  
  for (cc in 1:parms$chains) {
    
    # To store multiply imputed datasets (per chain)
    imp_DURR_dat <- vector("list", parms$iters)
      names(imp_DURR_dat) <- seq(1:parms$iters)
      
    # Storing imputate values for each iteration (per chain)
    imps <- lapply(1:p_imp, function(x) {
      matrix(data = NA, nrow = parms$iters, ncol = nr[x],
             dimnames = list(NULL, rownames(Z[!O[, x],]) ))
    })
    
    Zm <- init_dt_i(Z, missing_type(Z)) # reinitialize data
      for (i in 1:p_imp) imps[[i]][1, ] <- Zm[!O[, i], i]

    for (m in 2:parms$iters) {
      print(paste0("DURR - Chain: ", cc, "/", parms$chains, 
                   "; Iter: ", m, "/", parms$iters))
      for (j in 1:p_imp) {
        # 1. Generate bootstrap data set (for each j!)
        # Generate bootstrap sample
        idx_bs   <- sample(1:nrow(Zm), nrow(Zm), replace = TRUE)
        Zm_bs <- Zm[idx_bs, ]
        
        # Select data
        y_obs_bs <- z_j_obs  <- as.vector(Zm_bs[O[idx_bs, j], j])
        y_mis_bs <- zm_mj    <- as.vector(Zm_bs[!O[idx_bs, j], j]) # useless
        X_obs_bs <- Wm_j_obs <- as.matrix(Zm_bs[O[idx_bs, j], -j])
        X_mis_bs <- Wm_mj    <- as.matrix(Zm_bs[!O[idx_bs, j], -j])
        X_mis <- Wm_mj    <- as.matrix(Zm[!O[, j], -j])
        
        # J <- which(vmis_ind)[j]
        glmfam <- detect_family(Zm[, j])
        
        # 2. Fit regularized regression on bootstraped observed data
        ## Lasso
        if(reg_type == "lasso"){
          regu.mod <- rr_est_lasso(X = X_obs_bs, y = y_obs_bs, parms = parms)
        }
        
        ## Elastic net
        if(reg_type == "el"){
          regu.mod <- rr_est_elanet(X = X_obs_bs, y = y_obs_bs, parms = parms)
        }
        
        # 3. Predict zm_j (i.e. obtain imputations (imps))
        if(glmfam == "gaussian"){
          zm_j <- imp_gaus_DURR(regu.mod, X_obs_bs, y_obs_bs, X_mis, parms)
        }
        if(glmfam == "binomial"){
          zm_j <- imp_dich_DURR(regu.mod, X_obs_bs, y_obs_bs, X_mis, parms)
        }
        if(glmfam == "multinomial"){
          zm_j <- imp_multi_DURR(regu.mod, X_obs_bs, y_obs_bs, X_mis, parms)
        }
        
        # Append imputation
        Zm[is.na(Z[, j]), j] <- zm_j
        imps[[j]][m, ] <- zm_j # save iteration imputation
      }
      imp_DURR_dat[[m]] <- Zm
    }  
    imp_DURR_val[[cc]] <- imps
  }
  
  end.time <- Sys.time()
  
  return(list(dats = imp_DURR_dat[parms$keep_dt],
              imps = imp_DURR_val,
              time = difftime(end.time, 
                              start.time, 
                              units = "mins"))
  )
}
