### Title:    Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-05-19
### Notes:    function to impute data according to IURR method following 
###           reference papers Zhao Long 2016 (for univariate miss) and 
###           Deng et al 2016 (for multivariate miss).

impute_IURR <- function(Z, O, cond, reg_type="lasso", parms, perform = TRUE){
  ## For internals: run subroutines.R until impute_IURR, then
  # Z = Xy_mis
  # Z = Xy_input[, CIDX_all]
  # O = as.data.frame(!is.na(Z))            # matrix index of observed values
  # reg_type = c("el", "lasso")[2] # imputation model penality type
  
  ## Body
  if(perform == TRUE){
    
    tryCatch({
      
      p  <- ncol(Z) # number of variables [indexed with j]
      
      p_imp    <- sum(colMeans(O) < 1)
      p_imp_id <- names(which(colMeans(O) < 1))
      nr       <- colSums(!O[, colMeans(O) < 1])
      
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
        imps <- lapply(p_imp_id, function(x) {
          matrix(data = NA, nrow = parms$iters, ncol = nr[x],
                 dimnames = list(NULL, rownames(Z[!O[, x],]) ))
        })
        
        Zm <- init_dt_i(Z, missing_type(Z)) # reinitialize data
        imp_IURR_dat$`1` <- Zm
        for (i in 1:p_imp) imps[[i]][1, ] <- Zm[!O[, p_imp_id[i]], 
                                                p_imp_id[i]]
        
        for (m in 2:parms$iters) {
          print(paste0("IURR - Chain: ", cc, "/", parms$chains, 
                       "; Iter: ", m, "/", parms$iters))
          pb <- txtProgressBar(min = 0, max = p_imp, style = 3)
          
          for (j in 1:p_imp) {
            J <- which(colnames(Zm) %in% p_imp_id[j])
            # Select data
            y_obs <- z_j_obs  <- Zm[O[, J] == TRUE,  J]
            y_mis <- zm_mj    <- Zm[O[, J] == FALSE, J] # useless
            
            X_obs <- Zm[O[, J], -J] # Wm_j_obs
            X_obs <- model.matrix(~ ., X_obs)[, -1]
            X_mis <- Zm[!O[, J], -J] # Wm_mj
            X_mis <- model.matrix(~ ., X_mis)[, -1]
            
            glmfam <- detect_family(Zm[, J])
            
            # Fit regularized regression
            # Requirement for procedure is that selected n > p as the ML apporaches
            # used later in the estiamtion require this.
            IURR_fwd_cond <- FALSE
            while (IURR_fwd_cond == FALSE) {
              
              # Lasso
              if(reg_type == "lasso"){
                regu.mod <- rr_est_lasso(X = X_obs, y = y_obs, 
                                         parms = parms, fam = glmfam)
              }
              
              # Elastic net
              if(reg_type == "el"){
                regu.mod <- rr_est_elanet(X = X_obs, y = y_obs, 
                                          parms = parms, fam = glmfam)
              }
              
              IURR_fwd_cond <- length(y_obs) > sum(coef(regu.mod) != 0)
              success <- c(success, IURR_fwd_cond)
              
            }
            
            # 3. Predict zm_j (i.e. obtain imputations (imps))
            if(glmfam == "gaussian"){
              zm_j <- imp_gaus_IURR(model = regu.mod,
                                    X_tr  = X_obs, y_tr = y_obs,
                                    X_te  = X_mis, y_te = y_mis, 
                                    parms = parms)
            }
            if(glmfam == "binomial"){
              zm_j <- imp_dich_IURR(model = regu.mod,
                                    X_tr = X_obs, y_tr = y_obs,
                                    X_te = X_mis, 
                                    parms = parms)
            }
            if(glmfam == "multinomial"){
              zm_j <- imp_multi_IURR(regu.mod, X_obs_bs, y_obs_bs, X_mis, parms)
            }
            
            # Append imputation
            Zm[!O[, J], J] <- zm_j # update data
            imps[[j]][m, ] <- zm_j # save iteration imputation
            
            # Monitor Progress
            setTxtProgressBar(pb, j)
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
      
      
    }, error = function(e){
      err <- paste0("Original Error: ", e)
      print(err)
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
