### Title:    Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-06-09
### Notes:    Differs from version 1 in the code folder in what imputations
###           are stored: here we store the last imputed datasets of distinct
###           runs of the same number of iterations (similar to mice package)
###           while in the original, the last M imputed data sets after 
###           appropriate thinning are chosen for subsequent standard 
###           complete-data analysis (following DengEtAl 2016, see appendix).

impute_DURR_v2 <- function(Z, O, cond, reg_type="lasso", parms, perform = TRUE){
  ## Description
  # Z = Xy_mis
  # O <- as.data.frame(!is.na(Xy_mis))            # matrix index of observed values
  # reg_type = c("el", "lasso")[2] # imputation model penality type
  
  ## Body
  if(perform == TRUE){
    p  <- ncol(Z) # number of variables [indexed with j]
    
    p_imp <- length(parms$z_m_id) # variables with missing values
    r     <- colSums(O[parms$z_m_id]) # variable-wise count of observed values
    nr    <- colSums(!O[parms$z_m_id]) # variable-wise count of observed values
    
    # To store imputed values and check convergence
    imp_DURR_val <- vector("list", parms$chains)
     names(imp_DURR_val) <- seq(1:parms$chains)
      
    # To store multiply imputed datasets (per chain)
    imp_DURR_dat <- vector("list", parms$chains)
      names(imp_DURR_dat) <- seq(1:parms$chains)
      
    # Time performance
    start.time <- Sys.time()
    
    for (cc in 1:parms$chains) {
      # Storing imputated values for each iteration (per chain)
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
          y_obs_bs <- Zm_bs[O[idx_bs, j], j]  # z_j_obs
          y_mis_bs <- Zm_bs[!O[idx_bs, j], j] # zm_mj [useless]
          X_obs_bs <- as.matrix(Zm_bs[O[idx_bs, j], -j]) # Wm_j_obs
            X_obs_bs <- apply(X_obs_bs, 2, as.numeric) # makes dicho numbers
          X_mis <- Wm_mj    <- as.matrix(Zm[!O[, j], -j]) # Wm_mj
            X_mis <- apply(X_mis, 2, as.numeric)
          
          glmfam <- detect_family(Zm[, j])
          reg_type <- "lasso"
          
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
          Zm[is.na(Z[, j]), j] <- zm_j
          imps[[j]][m, ] <- zm_j # save iteration imputation thing
        }
      }
      colnames(Zm) <- colnames(Z)
      imp_DURR_dat[[cc]] <- Zm
      imp_DURR_val[[cc]] <- imps
    }
    
    end.time <- Sys.time()
    
    return(list(dats = imp_DURR_dat,
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
