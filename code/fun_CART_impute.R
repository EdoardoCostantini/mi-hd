### Title:    Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-06-26
### Notes:    function to impute data according to CART method
###           used in mice. This function differs from mice.impute.cart
###           becuase it allows to get different imputed datasets out
###           of 1 chain of imputations.

impute_CART <- function(Z, O, cond, parms, perform = TRUE){
  ## For internals: run subroutines.R until impute_CART, then
  # Z = Xy_mis
  # O <- as.data.frame(!is.na(Z))
  
  ## Body
  if(perform == TRUE){
    p  <- ncol(Z) # number of variables [indexed with ]
    nr       <- colSums(!O[, colMeans(O) < 1])

    # Craete a copy of Z to be processed
    Zm <- Z

    # To store imputed values and check convergence
    imp_CART_val <- vector("list", parms$chains)
    names(imp_CART_val) <- seq(1:parms$chains)
    
    # Time performance
    start.time <- Sys.time()
    
    for (cc in 1:parms$chains) {
      
      # To store multiply imputed datasets (per chain)
      imp_CART_dat <- vector("list", parms$iters)
      names(imp_CART_dat) <- seq(1:parms$iters)
      
      # Storing imputate values for each iteration (per chain)
      imps <- lapply(parms$z_m_id, function(x) {
        matrix(data = NA, nrow = parms$iters, ncol = nr[x],
               dimnames = list(NULL, rownames(Z[!O[, x],]) ))
      })
      
      Zm <- init_dt_i(Z, missing_type(Z)) # reinitialize data
      imp_CART_dat$`1` <- Zm
      for (i in 1:parms$zm_n) imps[[i]][1, ] <- Zm[!O[, parms$z_m_id[i]], 
                                              parms$z_m_id[i]]
      
      for (m in 2:parms$iters) {
        print(paste0("CART - Chain: ", cc, "/", parms$chains, 
                     "; Iter: ", m, "/", parms$iters))
        for (j in 1:parms$zm_n) {
          J <- which(colnames(Zm) %in% parms$z_m_id[j])
          # Select data
          y_obs <- z_j_obs  <- Zm[O[, J] == TRUE, J]
          y_mis <- zm_mj    <- Zm[O[, J] == FALSE, J] # useless
          X_obs <- Wm_j_obs <- as.matrix(Zm[O[, J] == TRUE, -J])
            X_obs <- apply(X_obs, 2, as.numeric) # makes dicho numbers
          X_mis <- Wm_mj    <- as.matrix(Zm[O[, J] == FALSE, -J])
            X_mis <- apply(X_mis, 2, as.numeric) # makes dicho numbers
          
          # Fit tree
          fit <- rpart::rpart(y_obs ~ ., 
                              data = as.data.frame(cbind(y_obs, X_obs)), 
                              method = "anova", 
                              control = rpart::rpart.control(minbucket = round(20/3), 
                                                                               cp = .01))
          # Generate imputations
          leafnr <- floor(as.numeric(row.names(fit$frame[fit$where, 
                                                         ])))
          fit$frame$yval <- as.numeric(row.names(fit$frame))
          nodes <- predict(object = fit, newdata = as.data.frame(X_mis))
          donor <- lapply(nodes, function(s) y_obs[leafnr == s])
          zm_j <- vapply(seq_along(donor), 
                         function(s) sample(donor[[s]], 
                                            1), 
                         numeric(1))
          
          # Append imputation
          Zm[!O[, J], J] <- zm_j # update data
          imps[[j]][m, ] <- zm_j # save iteration imputation
        }
        imp_CART_dat[[m]] <- Zm
      }
      imp_CART_val[[cc]] <- imps
    }
    
    end.time <- Sys.time()
    
    return(list(dats = imp_CART_dat[parms$keep_dt], # automatically form the last chain
                imps = imp_CART_val, # from all chains
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
