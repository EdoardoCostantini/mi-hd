# Title:    Imputation model using custom predictors
# Project:  Imputing High Dimensional Data
# Author:   Edoardo Costantini
# Created:  2023-03-19
# Modified: 2023-03-19

impute_stepwise <- function(Z, O, cond, direction = "forw", parms, perform = TRUE){
  ## Description
  # Z = Xy_mis
  # O = as.data.frame(!is.na(Z))            # matrix index of observed values
  # direction = c("forw", "back", "both")[1]
  
  ## Body
  if(perform == TRUE){
    
    tryCatch({
      
      p  <- ncol(Z) # number of variables [indexed with j]
      
      p_imp    <- sum(colMeans(O) < 1)
      p_imp_id <- names(which(colMeans(O) < 1))
      nr       <- colSums(!O[, colMeans(O) < 1])
      
      # To store imputed values and check convergence
      imp_stepwise_val <- vector("list", parms$chains)
        names(imp_stepwise_val) <- seq(1:parms$chains)
      
      # Time performance
      start.time <- Sys.time()
      
      for (cc in 1:parms$chains) {
        
        # To store multiply imputed datasets (per chain)
        imp_stepwise_dat <- vector("list", parms$iters)
          names(imp_stepwise_dat) <- seq(1:parms$iters)
        
        imps <- lapply(p_imp_id, function(x) {
          matrix(data = NA, nrow = parms$iters, ncol = nr[x],
                 dimnames = list(NULL, rownames(Z[!O[, x],]) ))
        })
        
        Zm <- init_dt_i(Z, missing_type(Z)) # reinitialize data
          imp_stepwise_dat$`1` <- Zm
          for (i in 1:p_imp) imps[[i]][1, ] <- Zm[!O[, p_imp_id[i]], 
                                                  p_imp_id[i]]

        for (m in 2:parms$iters) {
          # Monitor progress
          print(paste0("stepwise - Chain: ", cc, "/", parms$chains, 
                       "; Iter: ", m, "/", parms$iters, " at ", Sys.time()))
          
          for (j in 1:p_imp) {
            # Define index
            J <- which(colnames(Zm) %in% p_imp_id[j])

            # 1. Generate bootstrap data set (for each j!)
            idx_bs   <- sample(1:nrow(Zm), nrow(Zm), replace = TRUE)
            Zm_bs <- Zm[idx_bs, ]
            
            # Select data
            y_obs_bs <- Zm_bs[O[idx_bs, J], J]  # z_j_obs
            y_mis_bs <- Zm_bs[!O[idx_bs, J], J] # zm_mj [useless]
            
            X_obs_bs <- Zm_bs[O[idx_bs, J], -J] # Wm_j_obs
              X_obs_bs <- model.matrix(~., X_obs_bs)[, -1]
            X_mis <- Zm[!O[, J], -J] # Wm_mj
              X_mis <- model.matrix(~., X_mis)[, -1]

            # 2. Use stepwise to estimate model
            # Step For
            if (direction == "forw") {
              mod <- imputeR::stepForR(x = X_obs_bs, y = y_obs_bs)
            }
            if (direction == "back") {
              mod <- imputeR::stepBackR(x = X_obs_bs, y = y_obs_bs)
            }
            if (direction == "both") {
              mod <- imputeR::stepBothR(x = X_obs_bs, y = y_obs_bs)
            }

            # 3. Predict zm_j based on original data (not bootstrap)
            # (i.e. obtain imputations (imps))
            df <- max(1, mod$df.residual)
            s2hat <- sum(resid(mod)^2) / df
            
            # required to make sure that names are treated as in step function
            impdata <- data.frame(cbind(NA, X_mis))
            
            # Get predictions
            yhat <- as.matrix(cbind(1, impdata[, names(coef(mod))[-1]])) %*% coef(mod)
            
            # Make them imputations
            zm_j <- yhat + rnorm(nrow(X_mis), mean = 0, sd = sqrt(s2hat))

            # Append imputation
            Zm[is.na(Z[, J]), J] <- zm_j
            imps[[j]][m, ] <- zm_j # save iteration imputation thing

            # Monitor Progress
            # setTxtProgressBar(pb, j)
        }
          # close(pb)
          imp_stepwise_dat[[m]] <- Zm
        }  
        imp_stepwise_val[[cc]] <- imps
      }
      
      end.time <- Sys.time()
      
      return(list(dats = imp_stepwise_dat[parms$keep_dt],
                  imps = imp_stepwise_val,
                  time = difftime(end.time, 
                                  start.time, 
                                  units = "mins"))
      )
      
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
                time = NULL))
  }
}

impute_MICE_stepFor <- function(Z, perform = TRUE, parms = parms,
                           ridge = 0, eps = 0, threshold = 1){
  
  ## Input: 
  # @Z: dataset w/ missing values
  # @parms: the initialization object parms
  
  ## Example inputs
  # Z = Xy_mis # or
  # preds = parms$z_m_id
  # cond = conds[2,]
  # parms = parms
  
  ## output: 
  # - a list of chains imputed datasets at iteration iters
  # - per variable list of imputed values
  # - imputation run time
  
  ## body:
  if(perform == TRUE){
    
    tryCatch({
      
    start.time <- Sys.time()
    
    # Define predictor matrix based on quickpred
    pMat <- make.predictorMatrix(Z)

    # Define general imputation method
    methods <- make.method(
      Z,
      defaultMethod = c("imputeR.lmFun", "logreg", "polyreg", "polr")
    )

    # Define specifics of imputeR methods
    Fun <- list(
      z1 = imputeR::stepForR,
      z2 = imputeR::stepForR,
      z3 = imputeR::stepForR,
      z6 = imputeR::stepForR,
      z7 = imputeR::stepForR,
      z8 = imputeR::stepForR
    )

    # Impute
    imp_MIstepFor_mids <- mice(Z,
      predictorMatrix = pMat,
      m = parms$mice_ndt,
      maxit = parms$mice_iters,
      printFlag = TRUE,
      ridge = ridge,
      eps = eps,
      threshold = threshold,
      method = methods,
      Fun = Fun
    )

    # End time
    end.time <- Sys.time()

    # Store results
    imp_MIstepFor_dats <- mice::complete(imp_MIstepFor_mids, "all")
    imp_MIstepFor_imps <- imp_MIstepFor_mids$imp
    imp_MIstepFor_time <- difftime(end.time, start.time, units = "mins")
    
    return(list(dats = imp_MIstepFor_dats,
                imps = imp_MIstepFor_imps,
                time = imp_MIstepFor_time,
                mids = imp_MIstepFor_mids))
    
    
    }, error = function(e){
      err <- paste0("Original Error: ", e)
      print(err)
      return(list(dats = NULL,
                  imps = NULL,
                  time = NULL,
                  mids = NULL)
      )
    }
    )

  } else {
    return(list(dats = NULL,
                imps = NULL,
                time = NULL,
                mids = NULL))
  }
}
