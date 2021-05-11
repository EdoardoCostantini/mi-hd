### Title:    Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-05-19
### Notes:    function to impute data according to IURR method following 
###           reference papers Zhao Long 2016 (for univariate miss) and 
###           Deng et al 2016 (for multivariate miss).

impute_IURR <- function(Z, O, cond, reg_type = "lasso", parms, perform = TRUE){
  ## For internals: run subroutines.R until impute_IURR, then
  # Z = Xy_mis
  # Z = Xy_input[, CIDX_all]
  # O = as.data.frame(!is.na(Z))   # matrix index of observed values
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
      
      success <- 1 # stores TRUE/FALSE for success/failure of rr_est_function
      
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
        MLE_n_par <- vector("list", parms$iters)
        
        # Active Set Length Store
        store_AS_len <- matrix(NA, 
                               nrow = p_imp, 
                               ncol = parms$iters)
        
        for (m in 2:parms$iters) {
          print(paste0("IURR - Chain: ", cc, "/", parms$chains, 
                       "; Iter: ", m, "/", parms$iters, " at ", Sys.time()))
          
          # Progress for iteration m
          #pb <- txtProgressBar(min = 0, max = p_imp, style = 3)
          # Times
          for (j in 1:p_imp) {
            
            J <- which(colnames(Zm) %in% p_imp_id[j])
            # Select data
            y_obs <- Zm[O[, J] == TRUE,  J]
            y_mis <- Zm[O[, J] == FALSE, J] # useless
            
            X_obs <- as.matrix(Zm[O[, J], -J]) # Wm_j_obs
            # X_obs <- model.matrix(~ ., X_obs)[, -1]
            X_mis <- as.matrix(Zm[!O[, J], -J]) # Wm_mj
            # X_mis <- model.matrix(~ ., X_mis)[, -1]
            
            glmfam <- detect_family(Zm[, J])
            
            # 1. Define Active Set (by lasso)
            # Requirement for procedure is that selected n > p as the ML apporaches
            # used later in the estiamtion require this.
            cv_lasso <- cv.glmnet(x = X_obs, y = y_obs,
                                  family = glmfam,
                                  nfolds = 10, alpha = 1)
            
            # Extract Regularized Regression Coefs
            rr_coef <- as.matrix(coef(cv_lasso, s = "lambda.min"))[, 1]
            rr_coef_no0 <- names(rr_coef)[rr_coef != 0]
            rr_coef_noInt <- rr_coef_no0[-1]
            
            # Check if n-2 < p
            # -1 to estiamte sigma, -1 to esimate the intercept
            if((length(y_obs)-2) < length(rr_coef_noInt)){
              # Drop the smallest coefficients that make the system defined
              coef_sort <- sort(rr_coef[names(rr_coef) %in% rr_coef_noInt])
              coef_drop <- names(coef_sort[-(1:(length(y_obs)-2))])
              coef_drop_index <- c(1, which(rr_coef_no0 %in% coef_drop))
            } else {
              coef_drop_index <- 1 # just intercept
            }
            
            # Define Active Set
            AS <- rr_coef_no0[-coef_drop_index]
            store_AS_len[j, m] <- length(AS)
            
            # 2. Predict zm_j (i.e. obtain imputations (imps))
            # Define starting values
            if(identical(rr_coef_no0, "(Intercept)")){
              lm_fit <- lm(y_obs ~ 1)
              X_mle <- model.matrix(y_obs ~ 1)
            } else {
              lm_fit <- lm(y_obs ~ X_obs[, AS])
              X_mle  <- model.matrix(y_obs ~ X_obs[, AS])
              colnames(X_mle) <- str_replace(colnames(X_mle), ".*]+", "")
              b.estimated <- coef(lm_fit)[!is.na(coef(lm_fit))]
              b.names <- str_replace(names(b.estimated), ".*]+", "")
              X_mle  <- X_mle[, colnames(X_mle) %in% b.names]
              # Fix NAs when coefficinet cannot be estiamted because variable
              # is near constant
            }
            startV <- c(coef(lm_fit), sigma(lm_fit))

            # MLE estiamtes by Optimize loss function
            MLE_fit <- optim(startV,
                             .lm_loss,
                             method = "BFGS",
                             hessian = T,
                             y = y_obs, X = X_mle)
            theta <- MLE_fit$par
            OI <- solve(MLE_fit$hessian) # parameters cov matrix

            # Sample parameters for posterior predictive distribution
            pdraws_par <- MASS::mvrnorm(1, 
                                        mu = MLE_fit$par, 
                                        Sigma = OI)
            
            # 3. Predict zm_j (Sample posterior predictive distribution)
            if(identical(rr_coef_no0, "(Intercept)")){
              y_imp <- rnorm(n = nrow(X_mis),
                             mean = pdraws_par[1],
                             sd = pdraws_par[2])
            } else {
              X_ppd <- model.matrix( ~ X_mis[, AS]) # X for posterior pred dist
              colnames(X_ppd) <- str_replace(colnames(X_ppd), ".*]+", "")
              X_ppd  <- X_ppd[, colnames(X_ppd) %in% b.names]
              b_ppd <- pdraws_par[-length(pdraws_par)] # betas for posterior pred dist
              sigma_ppd <- tail(pdraws_par, 1) # sigma for posterior pred dist
              y_imp <- rnorm(n = nrow(X_mis),
                             mean = X_ppd %*% b_ppd,
                             sd = sigma_ppd)
            }
            # Append imputation
            Zm[!O[, J], J] <- y_imp # update data
            imps[[j]][m, ] <- y_imp # save iteration imputation
            # Monitor Progress
            #setTxtProgressBar(pb, j)
          }
          #close(pb)
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
                  # succ_ratio = mean(success))
                  succ_ratio = colMeans(store_AS_len, na.rm = T))
             
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
