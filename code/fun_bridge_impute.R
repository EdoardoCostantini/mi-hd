### Title:    Bayesian Ridge imputation
### Project:  Imputing High Dimensional Data
### Author:   Anonymized for peer review
### Created:  2020-05-19
### Modified: 2022-02-22

impute_BRIDGE <- function(Z, O, ridge_p, parms, perform = TRUE, robust = TRUE){
  
  # Prep data ---------------------------------------------------------------
  # robust = TRUE
  # cc <- 1
  # m <- 2
  # j <- 1

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
          print(paste0("bridge - Chain: ", cc, "/", parms$chains, 
                       "; Iter: ", m, "/", parms$iters, " at ", Sys.time()))
          # Loop across variables (cycle)
          for (j in 1:p_imp) {
            J <- which(colnames(Zm) %in% p_imp_id[j])
            rj <- !is.na(Z[, J])

            # Define observed ys
            zj_obs <- Zm[rj, J]

            # Separate predictor from dependent variable
            Zx <- Zm[, -J]

            # Create model matrix (no inercept)
            Zx <- model.matrix(~ ., Zx)[, -1]

            # Clean data
            if(robust == TRUE){
              # Find Dummies that have 99% of objservations in 1 category
              tabular <- apply(Zx, 2, table)

              # Select only dummies
              tabular <- tabular[sapply(tabular, length) == 2]

              # Vector of dummy names to discard
              dum.disc <- lapply(tabular, function(x) {
                x[1] / sum(x) > .95 | x[1] / sum(x) < 1-.95
              })
              dum.disc <- names(which(dum.disc == TRUE))
              Zx    <- Zx[, !colnames(Zx) %in% dum.disc]

              # Find collinear variables
              coll.vars <- find.collinear(Zx)
              Zx  <- Zx[, !colnames(Zx) %in% coll.vars]
            }

            # Center predictors (to avoid penalization of intercept)
            Zx <- apply(Zx, 2, function (col) col - mean(col))

            # Separate mis and obs parts of the data
            Zx_obs <- Zx[rj, ] # X values for missing ys
            Zx_mis <- Zx[!rj, ] # X values for missing ys

            # Add contants for intercept
            Zx_obs <- cbind(1, Zx_obs)
            Zx_mis <- cbind(1, Zx_mis)

            # Estimate values to define posterior distributions
            df <- max(length(zj_obs) - ncol(Zx_obs), 1)
            xtx <- crossprod(Zx_obs)        # get cross-product matrix
            pen <- ridge_p * diag(xtx)      # define vector of penalties
              pen[1] <- 0                   # make sure interecept is not penalized
            v <- solve(xtx + diag(pen))     # obtain penalised solution
            c <- t(zj_obs) %*% Zx_obs %*% v # obtain regression weights
            r <- zj_obs - Zx_obs %*% t(c)   # define residuals

            # Obtain posterior draws
            sigma_star <- sqrt(sum((r)^2) / rchisq(1, df))
            beta_star <- t(c) + (t(chol(v)) %*% rnorm(ncol(Zx_obs))) * sigma_star
            pdraw <- list(coef = c, beta = beta_star, sigma = sigma_star)

            # Obtain posterior predictive draws
            pdraw_zj_imp <- Zx_mis %*% pdraw$beta + rnorm(sum(!rj)) * pdraw$sigma

            # Append imputation (for next iteration)
            Zm[!rj, J] <- pdraw_zj_imp
            
            # Store imputations
            Imp.out[[j]][m, ]  <- pdraw_zj_imp
          }
          # only data from last chain will actually be saved
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

