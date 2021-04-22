### Title:    Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-05-19

impute_BLAS_hans <- function(Z, parms, perform = TRUE){
  
  # Prep data ---------------------------------------------------------------
  # Z = Xy_mis
  # Z = Xy_input[, CIDX_all]
  # O = as.data.frame(!is.na(Z))            # matrix index of observed values
  
  if(perform == TRUE){
    
    tryCatch({
      
      # Save NA action state
      default.na.action <- options('na.action')
      
      # Change to pass for model.matrix generation
      options(na.action = 'na.pass')
      
      # Generate model matrix
      Z_mm <- model.matrix(~., Z)[, -1]
      
      # Revert NA action
      options(na.action = default.na.action$na.action)
      
      # Get rid of constants (might happen for 1 unpopular dummy code)
      const <- names(which(apply(Z_mm, 2, var) == 0))
      Z_mm  <- Z_mm[, !colnames(Z_mm) %in% const]
      
      # Find Dummies that have 95% of observations in 1 category
      tabular <- apply(Z_mm, 2, table)
      
      # Select only dummies
      tabular <- tabular[sapply(tabular, length) == 2]
      
      # Vector of dummy names to discard
      dum.disc <- lapply(tabular, function(x) {
        x[1] / sum(x) > .95 | x[1] / sum(x) < 1-.95
      })
      dum.disc <- names(which(dum.disc == TRUE))
      Z_mm    <- Z_mm[, !colnames(Z_mm) %in% dum.disc]
      
      # Find collinear variables
      coll.vars <- find.collinear(Z_mm)
      Z_mm  <- data.frame(Z_mm[, !colnames(Z_mm) %in% coll.vars])
      
      # Generate observed/missing data locator object O
      O <- as.data.frame(!is.na(Z_mm))
      
      # Obtain data characteristics of model.matrix data version
      p        <- ncol(Z_mm)
      p_imp    <- sum(colMeans(O) < 1)
      p_imp_id <- names(which(colMeans(O) < 1))
      r        <- colSums(O[, colMeans(O) < 1])
      nr       <- colSums(!O[,colMeans(O) < 1])
      
      # To store imputed values and check convergence
      imp_blasso_val <- vector("list", parms$chains_bl)
        names(imp_blasso_val) <- seq(1:parms$chains_bl)
      
      # Time performance
      start.time <- Sys.time()
      
      # For one chain of imputatuions
      for (cc in 1:parms$chains_bl){
        Zm <- init_dt_i(Z_mm, missing_type(Z_mm)) # initialize data for each chain
        
        # To store multiply imputed datasets (in from the last chain)
        imp_blasso_dat <- vector("list", parms$iters_bl)
          names(imp_blasso_dat) <- seq(1:parms$iters_bl)
        imp_blasso_dat$`1` <- init_dt_i(Z, missing_type(Z)) # imputations from initialization
        
        # Empty storing objects for MCMC samples
        # For each, I need one slot per imputation model (target variable)
        
        # regression coefficients
        Beta.out <- lapply(1:p_imp, matrix, data = NA, nrow = parms$iters_bl,
                           ncol = ncol(Zm)-1)
        for (i in 1:p_imp) Beta.out[[i]][1, ] <- rep(0, ncol(Zm)-1)
        
        # Error variance
        Sig.out <- lapply(1:p_imp, matrix, data = NA, nrow = parms$iters_bl, 
                          ncol = 1)
        for (i in 1:p_imp) Sig.out[[i]][1, ] <- 1
        
        # tau parameter in Hans 2009 and 10 blasso model
        Tau.out <-  lapply(1:p_imp, matrix, data = NA, nrow = parms$iters_bl, 
                           ncol = 1)
        for (i in 1:p_imp) Tau.out[[i]][1, ] <- 1
        
        # phi parameter in Hans 2010 blasso model
        Phi.out <-  lapply(1:p_imp, matrix, data = NA, nrow = parms$iters_bl,
                           ncol = 1)
        for (i in 1:p_imp) Phi.out[[i]][1, ] <- .5
        
        # Imputed scores
        Imp.out <- lapply(p_imp_id, function(x) {
          matrix(data = NA, nrow = parms$iters_bl, ncol = nr[x],
                 dimnames = list(NULL, rownames(Zm[!O[, x],]) ))
        })
        for (i in 1:p_imp) Imp.out[[i]][1, ] <- Zm[!O[, p_imp_id[i]], 
                                                   p_imp_id[i]]
        
        # Latent variable for probit models
        lv.out <- lapply(1:p_imp, function(x) {
          matrix(NA, ncol = r[x], nrow = parms$iters_bl)
        })
        for (i in 1:p_imp) lv.out[[i]][1, ] <- 0
        
        # Loop across Iteration
        for (m in 2:parms$iters_bl) {
          Z_out <- Z # each iteration will output 1 dataset
          print(paste0("BLASSO - Chain: ",
                       cc, "/", parms$chains_bl, "; Iter: ",
                       m, "/", parms$iters_bl,
                       " at ", Sys.time()))
          # Loop across variables (cycle)
          for (j in 1:p_imp) {
            J <- which(colnames(Zm) %in% p_imp_id[j])
            zj_obs <- Zm[O[, J], J]
            zj_mis <- Zm[!O[, J], J]
            Z_obs  <- Zm[O[, J], -J]
            Z_mis  <- Zm[!O[, J], -J]
            
            # Scaled versions
            s_zj_obs <- scale(zj_obs)
            s_Z_obs  <- scale(Z_obs)
            s_Z_mis  <- scale(Z_mis,
                              center = colMeans(Z_obs),
                              scale = apply(Z_obs, 2, sd))
            
            # Storing draws  
            beta_m  <- Beta.out[[j]][m-1, ]
            sigma_m <- sqrt(Sig.out[[j]][m-1])
            tam_m   <- Tau.out[[j]][m-1]
            phi_m   <- Phi.out[[j]][m-1]
            
            # Detect type of vairbale for type of imputation
            fam <- detect_family(zj_obs)
            
            if(fam == "gaussian"){
              # step 1: update j-th imputation model parameters
              pdraw <- blasso::blasso.vs(Y = s_zj_obs, X = s_Z_obs,
                                         iters = 1,
                                         burn  = 0,
                                         beta  = beta_m, 
                                         beta.prior = "scaled", 
                                         sig2 = sigma_m, sig2prior = c(a = .1, b = .1),
                                         tau  = tam_m,   tauprior  = c(r = .01, s = .01), 
                                         phi  = phi_m,   phiprior  = c(h = 1, g = 1),
                                         fixsig = FALSE, fixtau = FALSE, fixphi = FALSE,
                                         noisy = FALSE)
              
              # step 2: sample from predictive distribution (predict on beta draws scale)
              s_pdraw_zj_imp <- rnorm(nrow(s_Z_mis),
                                      mean = (s_Z_mis %*% as.vector(pdraw$beta)), 
                                      sd = sqrt(pdraw$sig2) )
              # and rescale
              pdraw_zj_imp <- s_pdraw_zj_imp * sd(zj_obs) + mean(zj_obs)
            }
            if(fam == "binomial"){
              N1  <- sum(zj_obs == 1)     # Number of successes
              N0  <- length(zj_obs) - N1  # Number of failures
              
              # Draw latent variable z from its full conditional: z | \theta, y, X
              mu_lv <- Z_obs %*% beta_m # centers of truncated distributions for latent var
              
              lv.out[[j]][m, ][zj_obs == 0] <- rtruncnorm(N0, 
                                                          mean = mu_lv[zj_obs == 0], 
                                                          sd = sigma_m,
                                                          a = -Inf, b = 0)
              
              lv.out[[j]][m, ][zj_obs == 1] <- rtruncnorm(N1, 
                                                          mean = mu_lv[zj_obs == 1], 
                                                          sd = sigma_m,
                                                          a = 0, b = Inf)
              
              # Draw B_js
              pdraw <- blasso::blasso.vs(Y = lv.out[[j]][m, ], X = Z_obs,
                                         iters = 1,
                                         burn = 0,
                                         # Beta
                                         beta = beta_m, 
                                         beta.prior = "scaled", 
                                         # Sigma
                                         fixsig    = FALSE, # I should fix it to 1
                                         sig2      = sigma_m, 
                                         sig2prior = c(a = .1, b = .1),
                                         # Tau
                                         fixtau   = FALSE,
                                         tau      = tam_m,   
                                         tauprior = c(r = .01, s = .01), 
                                         # Phi
                                         fixphi   = FALSE,
                                         phi      = phi_m,
                                         phiprior = c(h = 1, g = 1),
                                         noisy = FALSE)
              # Sigma fix
              # 1. original blasso requires sigma prior for p > n
              # 2. original blasso samples tau last, using an updated sigma
              # 3. I modified blasso to sample sigma last. That way I can
              #    have 1 sample from blasso function where all values are
              #    sampled with sigma value = 1. I then want to overwrite
              #    the sampled sigma with value 1 so that next iteration
              #    will still use 1
              pdraw$sig2 <- 1 # I want this value to stay fixed at 1
              
              # Make prediction 
              P_y1 <- pnorm(Z_mis %*% as.vector(pdraw$beta))
              pdraw_zj_imp <- rbinom(nrow(Z_mis), 1, P_y1)
            }
            
            # Store results
            Beta.out[[j]][m, ] <- as.vector(pdraw$beta)
            Sig.out[[j]][m]    <- pdraw$sig2
            Tau.out[[j]][m]    <- ifelse(pdraw$tau > 0, pdraw$tau, Tau.out[[j]][m-1])
            # In unlikely scenario where Tau is smaller than 0, keep previous draw
            # [not sound]
            Phi.out[[j]][m]    <- pdraw$phi
            Imp.out[[j]][m, ]  <- pdraw_zj_imp
            
            # Append imputation (for next iteration)
            Zm[!O[, J], J] <- pdraw_zj_imp
          }
          # Replace imputed columns with imputated ones
          Z_out[, p_imp_id] <- Zm[, p_imp_id]
          imp_blasso_dat[[m]] <- Z_out
        }
        imp_blasso_val[[cc]] <- Imp.out
      }
      
      end.time <- Sys.time()
      
      return(list(dats = imp_blasso_dat[parms$keep_dt_bl],
                  imps = imp_blasso_val,
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
                time = NULL,
                succ_ratio = NULL)
    )
  }
}
