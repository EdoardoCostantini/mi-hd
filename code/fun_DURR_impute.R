### Title:    Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-05-19
### Notes:    function to impute data according to DURR method following 
###           reference papers Zhao Long 2016 (for univariate miss) and 
###           Deng et al 2016 (for multivariate miss)

impute_DURR <- function(Z, O, cond, reg_type = "lasso", parms, perform = TRUE){
  ## Description
  # Z = Xy_mis_DA_fill[, CIDX_all]
  # O <- as.data.frame(!is.na(Z))            # matrix index of observed values
  # reg_type = c("el", "lasso")[2] # imputation model penality type
  
  ## Body
  if(perform == TRUE){
    p  <- ncol(Z) # number of variables [indexed with j]
    nr <- colSums(!O[, colMeans(O) < 1])
    
    # Craete a copy of Z to be processed
    Zm   <- Z
    
    # To store imputed values and check convergence
    imp_DURR_val <- vector("list", parms$chains)
      names(imp_DURR_val) <- seq(1:parms$chains)
    
    # Time performance
    start.time <- Sys.time()
    
    for (cc in 1:parms$chains) {
      
      # To store multiply imputed datasets (per chain)
      imp_DURR_dat <- vector("list", parms$iters)
        names(imp_DURR_dat) <- seq(1:parms$iters)
      
      # Create storing object for imputated values
      imps <- lapply(parms$z_m_id, function(x) {
        matrix(data = NA, nrow = parms$iters, ncol = nr[x],
               dimnames = list(NULL, rownames(Z[!O[, x],]) ))
      })
      
      # Fill in and store initial imputations
      Zm[, parms$z_m_id] <- init_dt_i(Z[, parms$z_m_id], 
                                      missing_type(Z[, parms$z_m_id]))
        imp_DURR_dat$`1` <- Zm
      for (i in 1:parms$zm_n) imps[[i]][1, ] <- Zm[!O[, parms$z_m_id[i]], 
                                                   parms$z_m_id[i]]
        
      for (m in 2:parms$iters) {
        print(paste0("DURR - Chain: ", cc, "/", parms$chains, 
                     "; Iter: ", m, "/", parms$iters))
        for (j in 1:parms$zm_n) {
          J <- which(colnames(Zm) %in% parms$z_m_id[j])
          # 1. Generate bootstrap data set (for each j!)
          # Generate bootstrap sample
          idx_bs <- sample(1:nrow(Zm), nrow(Zm), replace = TRUE)
          Zm_bs  <- Zm[idx_bs, ]
          
          # Select data
          y_obs_bs <- Zm_bs[O[idx_bs, J], J]  # z_j_obs
          y_mis_bs <- Zm_bs[!O[idx_bs, J], J] # zm_mj [useless]
          X_obs_bs <- as.matrix(Zm_bs[O[idx_bs, J], -J]) # Wm_j_obs
            X_obs_bs <- apply(X_obs_bs, 2, as.numeric) # makes dicho numbers
          X_mis <- as.matrix(Zm[!O[, J], -J]) # Wm_mj
            X_mis <- apply(X_mis, 2, as.numeric)
          
          glmfam <- detect_family(Zm[, J])
          
          # 2. Fit regularized regression on bootstraped observed data
          ## Lasso
          if(reg_type == "lasso"){
            regu.mod <- rr_est_lasso(X = X_obs_bs, y = y_obs_bs, 
                                     parms = parms, fam = glmfam)
          }

          ## Elastic net
          if(reg_type == "el"){
            regu.mod <- rr_est_elanet(X = X_obs_bs, y = y_obs_bs, 
                                      parms = parms, fam = glmfam)
          }
          
          # 3. Predict zm_j based on onrignal data (not bootsrap)
          # (i.e. obtain imputations (imps))
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
          Zm[is.na(Z[, J]), J] <- zm_j
          imps[[j]][m, ] <- zm_j # save iteration imputation thing
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
  } else {
    return(list(dats = NULL,
                imps = NULL,
                time = NULL))
  }
}
