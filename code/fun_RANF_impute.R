### Title:    Imputing High Dimensional Data
### Author:   Anonymized for peer review
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
    
    tryCatch({
    
    p  <- ncol(Z) # number of variables [indexed with j]
    
    p_imp    <- sum(colMeans(O) < 1)
    p_imp_id <- names(which(colMeans(O) < 1))
    nr       <- colSums(!O[, colMeans(O) < 1])
    
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
      imps <- lapply(p_imp_id, function(x) {
        matrix(data = NA, nrow = parms$iters, ncol = nr[x],
               dimnames = list(NULL, rownames(Z[!O[, x],]) ))
      })
      
      Zm <- init_dt_i(Z, missing_type(Z)) # reinitialize data
      imp_RANF_dat$`1` <- Zm
      for (i in 1:p_imp) imps[[i]][1, ] <- Zm[!O[, p_imp_id[i]], 
                                              p_imp_id[i]]
      
      for (m in 2:parms$iters) {
        print(paste0("RANF - Chain: ", cc, "/", parms$chains, 
                     "; Iter: ", m, "/", parms$iters, " at ", Sys.time()))
        for (j in 1:p_imp) {
          J <- which(colnames(Zm) %in% p_imp_id[j])
          # Select data
          y_obs <- z_j_obs  <- Zm[O[, J] == TRUE, J]
          y_mis <- zm_mj    <- Zm[O[, J] == FALSE, J] # useless
          X_obs <- Zm[O[, J] == TRUE, -J]
          X_mis <- Zm[O[, J] == FALSE, -J]
            
          # Fit random forest
          suppressWarnings(
            # In exp4, where ordered items w/ 4 categories need to be imputed
            # having the item as numeric passed to the randomForest function
            # produces a warning asking whether you are sure you want to use
            # regression. To use a decision tree you would need to pass the
            # y_obs argument as a factor. However, for consistency we stick
            # to the numeric treatment and suppress the warnings that come
            # out of it
            forest <- sapply(seq_len(parms$rfntree), 
                             FUN = function(s) onetree(X_obs,
                                                       X_mis, 
                                                       y_obs))
          )
          # Generate imputations
          zm_j <- apply(forest, 
                        1, 
                        # Each observation with a missing value has its own pool
                        # of donors, which was found by fitting 10 different single
                        # trees. Possible donors from all different trees are
                        # used to make up a larger pool of donors from which to sample
                        # 1 imputed value
                        FUN = function(s) sample(unlist(s), 1))
          
          # Append imputation
          Zm[!O[, J], J] <- zm_j # update data
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
    
    ### END TRYCATCH EXPRESSION
    }, error = function(e){
      err <- paste0("Original Error: ", e)
      print(err)
      return(list(dats = NULL,
                  imps = NULL,
                  time = NULL)
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
