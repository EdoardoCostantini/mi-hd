### Title:    Bayesian Ridge imputation
### Project:  Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-05-19

impute_BRIDGE <- function(Z, O, ridge_p, parms, perform = TRUE){
  
  # Prep data ---------------------------------------------------------------
  # Z = Xy_mis
  # O = as.data.frame(!is.na(Xy_mis))            # matrix index of observed values
  # ridge_p = cond$ridge
  if(perform == TRUE){
    
    tryCatch({
      p  <- ncol(Z) # number of variables [indexed with J]
      
      p_imp    <- sum(colMeans(O) < 1)
      p_imp_id <- names(which(colMeans(O) < 1))
      nr       <- colSums(!O[, colMeans(O) < 1])
      
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
        Imp.out <- lapply(p_imp_id, function(x) {
          matrix(data = NA, nrow = parms$iters, ncol = nr[x],
                 dimnames = list(NULL, rownames(Zm[!O[, x],]) ))
        })
        for (i in 1:p_imp) Imp.out[[i]][1, ] <- Zm[!O[, p_imp_id[i]], 
                                                   p_imp_id[i]]
        
        # Loop across Iteration
        for (m in 2:parms$iters) {
          print(paste0("bridge - Chain: ", cc, "/", parms$chains, "; Iter: ", m, "/", parms$iters))
          # Loop across variables (cycle)
          for (j in 1:p_imp) {
            J <- which(colnames(Zm) %in% p_imp_id[j])
            zj_obs <- Zm[O[, J], J]
            zj_mis <- Zm[!O[, J], J]
            Z_obs <- as.matrix(Zm[O[, J], -J])
            Z_obs <- apply(Z_obs, 2, as.numeric) # makes dicho numbers
            Z_mis <- as.matrix(Zm[!O[, J], -J])
            Z_mis <- apply(Z_mis, 2, as.numeric) # makes dicho numbers
            
            # Obtain posterior draws for all paramters of interest
            pdraw <- .norm.draw(y       = Zm[, J], 
                                ry      = O[, J], 
                                x       = as.matrix(Zm[, -J]), 
                                ls.meth = "ridge", ridge = ridge_p)
            
            # Obtain posterior predictive draws
            pdraw_zj_imp <- Z_mis %*% pdraw$beta + rnorm(sum(!O[, J])) * pdraw$sigma
            
            # Append imputation (for next iteration)
            Zm[!O[, J], J] <- pdraw_zj_imp
            
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

