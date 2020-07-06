### Title:    Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-05-19

impute_BLAS_hans <- function(Z, O, parms, perform = TRUE){
  
  # Prep data ---------------------------------------------------------------
  # Z = Xy_mis
  # O = as.data.frame(!is.na(Xy_mis))            # matrix index of observed values
  if(perform == TRUE){
  p  <- ncol(Z) # number of variables [indexed with j]

  p_imp <- length(parms$z_m_id) # variables with missing values
  r  <- colSums(O[1:p_imp]) # variable-wise count of observed values
  nr  <- colSums(!O[1:p_imp]) # variable-wise count of observed values
  
  # To store imputed values and check convergence
  imp_blasso_val <- vector("list", parms$chains)
    names(imp_blasso_val) <- seq(1:parms$chains)

  # Time performance
  start.time <- Sys.time()
  
  # For one chain of imputatuions
  for (cc in 1:parms$chains) {
    
    # To store multiply imputed datasets (in from the last chain)
    imp_blasso_dat <- vector("list", parms$iters)
      names(imp_blasso_dat) <- seq(1:parms$iters)
    
    Zm <- init_dt_i(Z, missing_type(Z)) # initialize data for each chain
      imp_blasso_dat$`1` <- Zm
    # Empty storing objects for MCMC samples
    # For each, I need one slot per imputaiton model
      
    # regression coefficients
    Beta.out <- lapply(1:p_imp, matrix, data= NA, nrow=parms$iters, ncol=ncol(Z)-1)
    for (i in 1:p_imp) Beta.out[[i]][1, ] <- rep(0, ncol(Z)-1)
    
    # Error variance
    Sig.out <- lapply(1:p_imp, matrix, data = NA, nrow = parms$iters, ncol = 1)
    for (i in 1:p_imp) Sig.out[[i]][1, ] <- 1
    
    # tau parameter in Hans 2009 and 10 blasso model
    Tau.out <-  lapply(1:p_imp, matrix, data = NA, nrow = parms$iters, ncol = 1)
    for (i in 1:p_imp) Tau.out[[i]][1, ] <- 1
    
    # phi parameter in Hans 2010 blasso model
    Phi.out <-  lapply(1:p_imp, matrix, data = NA, nrow = parms$iters, ncol = 1)
    for (i in 1:p_imp) Phi.out[[i]][1, ] <- .5
    
    # Imputed scores
    Imp.out <- lapply(1:p_imp, function(x) {
      matrix(data = NA, nrow = parms$iters, ncol = nr[x],
             dimnames = list(NULL, rownames(Zm[!O[, x],]) ))
    })
    for (i in 1:p_imp) Imp.out[[i]][1, ] <- Zm[!O[, i], i]
    
    # Latent variable for probit models
    lv.out <- lapply(1:p_imp, function(x) {
      matrix(NA, ncol = r[x], nrow = parms$iters)
    })
    for (i in 1:p_imp) lv.out[[i]][1, ] <- 0
    
    # Loop across Iteration
    for (m in 2:parms$iters) {
      print(paste0("BLASSO - Chain: ", cc, "/", parms$chains, "; Iter: ", m, "/", parms$iters))
      # Loop across variables (cycle)
      for (j in 1:p_imp) {
        zj_obs <- Zm[O[, j], j]
        zj_mis <- Zm[!O[, j], j]
        Z_obs <- as.matrix(Zm[O[, j], -j])
          Z_obs <- apply(Z_obs, 2, as.numeric) # makes dicho numbers
        Z_mis <- as.matrix(Zm[!O[, j], -j])
          Z_mis <- apply(Z_mis, 2, as.numeric) # makes dicho numbers
          
        beta_m  <- Beta.out[[j]][m-1,] #rep(0, ncol(Z_obs)) # starting vbalues for gibbs sampler
        sigma_m <- sqrt(Sig.out[[j]][m-1])
        tam_m   <- Tau.out[[j]][m-1]
        phi_m   <- Phi.out[[j]][m-1]
        
        # Detect type of vairbale for type of imputation
        fam <- detect_family(zj_obs)
        
        if(fam == "gaussian"){
          # step 1: update j-th imputation model parameters
          pdraw <- blasso::blasso.vs(Y = zj_obs, X = Z_obs,
                                     iters =1,
                                     burn = 0,
                                     beta = beta_m, 
                                     beta.prior = "scaled", 
                                     sig2 = sigma_m, sig2prior = c(a = .1, b = .1),
                                     tau  = tam_m,   tauprior  = c(r = .01, s = .01), 
                                     phi  = phi_m,   phiprior  = c(h = 1, g = 1),
                                     fixsig = FALSE, fixtau = FALSE, fixphi = FALSE,
                                     noisy = FALSE)

          # step 2: update imputed values (sample from predictive distribution)
          pdraw_zj_imp <- rnorm(nrow(Z_mis),
                                mean = (Z_mis %*% as.vector(pdraw$beta)), 
                                sd = sqrt(pdraw$sig2) )
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
        Tau.out[[j]][m]    <- pdraw$tau
        Phi.out[[j]][m]    <- pdraw$phi
        Imp.out[[j]][m, ]  <- pdraw_zj_imp
        
        # Append imputation (for next iteration)
        Zm[!O[, j], j] <- pdraw_zj_imp
      }
      # only data from last chain will actaully be saved
      imp_blasso_dat[[m]] <- Zm 
    }
    imp_blasso_val[[cc]] <- Imp.out
  }

  end.time <- Sys.time()
  
  return(list(dats = imp_blasso_dat[parms$keep_dt],
              imps = imp_blasso_val,
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

