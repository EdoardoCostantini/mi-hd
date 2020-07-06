### Title:    Bayesian Ridge imputation
### Project:  Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-05-19

impute_BRIDGE <- function(Z, O, parms, perform = TRUE){
  
  # Prep data ---------------------------------------------------------------
  # Z = Xy_mis
  # O = as.data.frame(!is.na(Xy_mis))            # matrix index of observed values
  if(perform == TRUE){
    
    tryCatch({
      p  <- ncol(Z) # number of variables [indexed with j]
      
      p_imp <- length(parms$z_m_id) # variables with missing values
      r  <- colSums(O[1:p_imp]) # variable-wise count of observed values
      nr  <- colSums(!O[1:p_imp]) # variable-wise count of observed values
      
      # To store imputed values and check convergence
      imp_bridge_val <- vector("list", parms$chains)
      names(imp_bridge_val) <- seq(1:parms$chains)
      
      # Time performance
      start.time <- Sys.time()
      
      # For one chain of imputatuions
      for (cc in 1:parms$chains) {
        
        # To store multiply imputed datasets (in from the last chain)
        imp_bridge_dat <- vector("list", parms$iters)
        names(imp_bridge_dat) <- seq(1:parms$iters)
        
        Zm <- init_dt_i(Z, missing_type(Z)) # initialize data for each chain
        imp_bridge_dat$`1` <- Zm
        # Empty storing objects for MCMC samples
        
        # Imputed scores
        Imp.out <- lapply(1:p_imp, function(x) {
          matrix(data = NA, nrow = parms$iters, ncol = nr[x],
                 dimnames = list(NULL, rownames(Zm[!O[, x],]) ))
        })
        for (i in 1:p_imp) Imp.out[[i]][1, ] <- Zm[!O[, i], i]
        
        # Loop across Iteration
        for (m in 2:parms$iters) {
          print(paste0("bridge - Chain: ", cc, "/", parms$chains, "; Iter: ", m, "/", parms$iters))
          # Loop across variables (cycle)
          for (j in 1:p_imp) {
            zj_obs <- Zm[O[, j], j]
            zj_mis <- Zm[!O[, j], j]
            Z_obs <- as.matrix(Zm[O[, j], -j])
            Z_obs <- apply(Z_obs, 2, as.numeric) # makes dicho numbers
            Z_mis <- as.matrix(Zm[!O[, j], -j])
            Z_mis <- apply(Z_mis, 2, as.numeric) # makes dicho numbers
            
            # Obtain posterior draws for all paramters of interest
            pdraw <- .norm.draw(y       = Zm[, j], 
                                ry      = O[, j], 
                                x       = as.matrix(Zm[, -j]), 
                                ls.meth = "ridge", ridge = parms$ridge)
            
            # Obtain posterior predictive draws
            pdraw_zj_imp <- Z_mis %*% pdraw$beta + rnorm(sum(!O[, j])) * pdraw$sigma
            
            # Append imputation (for next iteration)
            Zm[!O[, j], j] <- pdraw_zj_imp
            
            # Store imputations
            Imp.out[[j]][m, ]  <- pdraw_zj_imp
          }
          # only data from last chain will actaully be saved
          imp_bridge_dat[[m]] <- Zm 
        }
        imp_bridge_val[[cc]] <- Imp.out
      }
      
      end.time <- Sys.time()
      
      return(list(dats = imp_bridge_dat[parms$keep_dt],
                  imps = imp_bridge_val,
                  time = difftime(end.time, 
                                  start.time, 
                                  units = "mins"))
      )
    }, error = function(e){
      return(list(dats = NULL,
                  imps = NULL,
                  time = NULL,
                  succ_ratio = NULL)
      )
    }
    )

  } else {
    return(list(dats = NULL,
                imps = NULL,
                time = NULL,
                succ_ratio = NULL)
    )
  }
}

