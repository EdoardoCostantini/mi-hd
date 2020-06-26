### Title:    Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-06-26
### Notes:    function to impute data according to Random Forest method
###           used in mice. This function differs from mice.impute.rf
###           becuase it allows to get different imputed datasets out
###           of 1 chain of imputations

impute_RANF <- function(Z, O, cond, parms, perform = TRUE){
  ## For internals: run subroutines.R until impute_RANF, then
  # Z = Xy_mis
  # O <- as.data.frame(!is.na(Xy_mis))            # matrix index of observed values
  
  ## Body
  if(perform == TRUE){
    p  <- ncol(Z) # number of variables [indexed with j]
    
    p_imp <- length(parms$z_m_id) # variables with missing values
    r     <- colSums(O[1:p_imp]) # variable-wise count of observed values
    nr    <- colSums(!O[1:p_imp]) # variable-wise count of observed values
    
    # To store imputed values and check convergence
    imp_RANF_val <- vector("list", parms$chains)
    names(imp_RANF_val) <- seq(1:parms$chains)
    
    # Time performance
    start.time <- Sys.time()
    
    for (cc in 1:parms$chains) {
      
      # To store multiply imputed datasets (per chain)
      imp_RANF_dat <- vector("list", parms$iters)
      names(imp_RANF_dat) <- seq(1:parms$iters)
      
      # Storing imputate values for each iteration (per chain)
      imps <- lapply(1:p_imp, function(x) {
        matrix(data = NA, nrow = parms$iters, ncol = nr[x],
               dimnames = list(NULL, rownames(Z[!O[, x],]) ))
      })
      
      Zm <- init_dt_i(Z, missing_type(Z)) # reinitialize data
      imp_RANF_dat$`1` <- Zm
      for (i in 1:p_imp) imps[[i]][1, ] <- Zm[!O[, i], i]
      
      for (m in 2:parms$iters) {
        print(paste0("RANF - Chain: ", cc, "/", parms$chains, 
                     "; Iter: ", m, "/", parms$iters))
        for (j in 1:p_imp) {
          # Select data
          y_obs <- z_j_obs  <- Zm[O[, j] == TRUE, j]
          y_mis <- zm_mj    <- Zm[O[, j] == FALSE, j] # useless
          X_obs <- Wm_j_obs <- as.matrix(Zm[O[, j] == TRUE, -j])
            X_obs <- apply(X_obs, 2, as.numeric) # makes dicho numbers
          X_mis <- Wm_mj    <- as.matrix(Zm[O[, j] == FALSE, -j])
            X_mis <- apply(X_mis, 2, as.numeric) # makes dicho numbers
          
          # Fit random forest
          forest <- sapply(seq_len(parms$rfntree), 
                           FUN = function(s) onetree(X_obs,
                                                     X_mis, 
                                                     y_obs))
          
          # Generate imputations
          zm_j <- apply(forest, 1, 
                        FUN = function(s) sample(unlist(s), 1))
          
          # Append imputation
          Zm[!O[, j], j] <- zm_j # update data
          imps[[j]][m, ] <- zm_j # save iteration imputation
        }
        imp_RANF_dat[[m]] <- Zm
      }
      imp_RANF_val[[cc]] <- imps
    }
    
    end.time <- Sys.time()
    
    return(list(dats = imp_RANF_dat[parms$keep_dt], # automatically form the last chain
                imps = imp_RANF_val, # from all chains
                time = difftime(end.time, 
                                start.time, 
                                units = "mins"))
    )
  } else {
    return(list(dats = NULL,
                imps = NULL,
                time = NULL,
                succ_ratio = NULL)
    )
  }
  
}
