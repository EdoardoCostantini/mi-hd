### Title:    Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-05-19
### Notes:    function to impute data according to IURR method following 
###           reference papers Zhao Long 2016 (for univariate miss) and 
###           Deng et al 2016 (for multivariate miss).

impute_IURR <- function(Z, O, cond, reg_type="lasso", parms, perform = TRUE){
  ## For internals: run subroutines.R until impute_IURR, then
  # Z = Xy_mis
  # O <- as.data.frame(!is.na(Xy_mis))            # matrix index of observed values
  # reg_type = c("el", "lasso")[2] # imputation model penality type
  
  ## Body
  if(perform == TRUE){
    p  <- ncol(Z) # number of variables [indexed with j]
    
    p_imp <- length(parms$z_m_id) # variables with missing values
    r  <- colSums(O[parms$z_m_id]) # variable-wise count of observed values
    nr  <- colSums(!O[parms$z_m_id]) # variable-wise count of observed values
    
    # To store imputed values and check convergence
    imp_IURR_val <- vector("list", parms$chains)
      names(imp_IURR_val) <- seq(1:parms$chains)
  
    # Time performance
    start.time <- Sys.time()
    
    success <- NULL # stores TRUE/FALSE for success/failure of rr_est_function
    
    for (cc in 1:parms$chains) {
    
      # To store multiply imputed datasets (per chain)
      imp_IURR_dat <- vector("list", parms$iters)
        names(imp_IURR_dat) <- seq(1:parms$iters)
        
      # Storing imputate values for each iteration (per chain)
      imps <- lapply(1:p_imp, function(x) {
        matrix(data = NA, nrow = parms$iters, ncol = nr[x],
               dimnames = list(NULL, rownames(Z[!O[, x],]) ))
      })
      
      Zm <- init_dt_i(Z, missing_type(Z)) # reinitialize data
        imp_IURR_dat$`1` <- Zm
        for (i in 1:p_imp) imps[[i]][1, ] <- Zm[!O[, i], i]
      
      for (m in 2:parms$iters) {
        print(paste0("IURR - Chain: ", cc, "/", parms$chains, 
                     "; Iter: ", m, "/", parms$iters))
        for (j in 1:p_imp) {
          # Select data
          y_obs <- z_j_obs  <- as.vector(Zm[O[, j] == TRUE, j])
          y_mis <- zm_mj    <- as.vector(Zm[O[, j] == FALSE, j]) # useless
          X_obs <- Wm_j_obs <- as.matrix(Zm[O[, j] == TRUE, -j])
          X_mis <- Wm_mj    <- as.matrix(Zm[O[, j] == FALSE, -j])
          
          # Fit regularized regression
          # Requirement for procedure is that selected n > p as the ML apporaches
          # used later in the estiamtion require this.
          IURR_fwd_cond <- FALSE
          while (IURR_fwd_cond == FALSE) {
            
            # Lasso
            if(reg_type == "lasso"){
              regu.mod <- rr_est_lasso(X = X_obs, y = y_obs, parms = parms)
            }
            
            # Elastic net
            if(reg_type == "el"){
              regu.mod <- rr_est_elanet(X = X_obs, y = y_obs, parms = parms)
            }
            
            IURR_fwd_cond <- length(z_j_obs) > sum(coef(regu.mod) != 0)
            success <- c(success, IURR_fwd_cond)
            
          }
          
          # Impute
          zm_j <- imp_gaus_IURR(model = regu.mod,
                                X_tr = Wm_j_obs, y_tr = z_j_obs,
                                X_te = Wm_mj, y_te = zm_mj, 
                                parms = parms)
          
          # Append imputation
          Zm[!O[, j], j] <- zm_j # update data
          imps[[j]][m, ] <- zm_j # save iteration imputation
        }
      imp_IURR_dat[[m]] <- Zm
      }
      imp_IURR_val[[cc]] <- imps
    }
    
    end.time <- Sys.time()
    
    return(list(dats = imp_IURR_dat[parms$keep_dt], # automatically form the last chain
                imps = imp_IURR_val, # from all chains
                time = difftime(end.time, 
                                start.time, 
                                units = "mins"),
                succ_ratio = mean(success))
    )
  } else {
    return(list(dats = NULL,
                imps = NULL,
                time = NULL,
                succ_ratio = NULL)
    )
  }
  
}
