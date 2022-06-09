### Title:    imputeHD-comp impute w/ Regularized Frequentiest Regressions DURR
### Author:   Anonymized for peer review
### Created:  2020-05-04
### Notes:    Function can be used only for univariate missingness imputation

impute_BLAS_hans <- function(Z, O, parms){
  
  # Prep data ---------------------------------------------------------------
  # Z = Xy_mis
  # O = as.data.frame(!is.na(Xy_mis))            # matrix index of observed values
  
  p  <- ncol(Z) # number of variables [indexed with j]

  p_imp <- length(parms$z_m_id) # variables with missing values
  r  <- colSums(O[parms$z_m_id]) # variable-wise count of observed values
  nr  <- colSums(!O[parms$z_m_id]) # variable-wise count of observed values
  
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
    # Empty storing objects for MCMC samples
    # For each, I need one slot per imputaiton model
    Beta.out <- lapply(1:p_imp, matrix, data= NA, nrow=parms$iters, ncol=ncol(Z)-1)
    for (i in 1:p_imp) Beta.out[[i]][1, ] <- rep(0, ncol(Z)-1)
    
    Sig.out <- lapply(1:p_imp, matrix, data = NA, nrow = parms$iters, ncol = 1)
    for (i in 1:p_imp) Sig.out[[i]][1, ] <- 1
    
    Tau.out <-  lapply(1:p_imp, matrix, data = NA, nrow = parms$iters, ncol = 1)
    for (i in 1:p_imp) Tau.out[[i]][1, ] <- 1
    
    Phi.out <-  lapply(1:p_imp, matrix, data = NA, nrow = parms$iters, ncol = 1)
    for (i in 1:p_imp) Phi.out[[i]][1, ] <- .5
    
    Imp.out <-  lapply(1:p_imp, function(x) {
      matrix(data = NA, nrow = parms$iters, ncol = nr[x],
             dimnames = list(NULL, rownames(Zm[!O[, x],]) ))
    })
    for (i in 1:p_imp) Imp.out[[i]][1, ] <- Zm[!O[, i], i]
    
    for (m in 2:parms$iters) {
      print(paste0("BLASSO - Chain: ", cc, "/", parms$chains, "; Iter: ", m, "/", parms$iters))
      for (j in 1:p_imp) {
        zj_obs <- as.vector(Zm[O[, j], j])
        zj_mis <- as.vector(Zm[!O[, j], j])
        Z_obs <- as.matrix(Zm[O[, j], -j])
        Z_mis <- as.matrix(Zm[!O[, j], -j])
        
        beta_m  <- Beta.out[[j]][m-1,] #rep(0, ncol(Z_obs)) # starting vbalues for gibbs sampler
        sigma_m <- sqrt(Sig.out[[j]][m-1])
        tam_m   <- Tau.out[[j]][m-1]
        phi_m   <- Phi.out[[j]][m-1]
        
        # step 1: update j-th imputation model parameters
        par_pdraws <- blasso::blasso.vs(Y = zj_obs, X = Z_obs,
                                        iters =1,
                                        burn = 0,
                                        beta = beta_m, 
                                        beta.prior = "scaled", 
                                        sig2 = sigma_m, sig2prior = c(a = .1, b = .1),
                                        tau  = tam_m,   tauprior  = c(r = .01, s = .01), 
                                        phi  = phi_m,   phiprior  = c(h = 1, g = 1),
                                        fixsig = FALSE, fixtau = FALSE, fixphi = FALSE,
                                        noisy = FALSE)
        # Gibbs sampler needs to nest the 1 sample per variable per iteration
        # check if blasso.vs akllows for taking 1 sample, then udpate everything 
        # absed on it. Imputation of every vairable is the target of the gibbs sampler
        # most inernal then oop over variable need to impte, that is one gibss sampler
        # multiple chians here are just for convergnece checks. Sample from the pool cof converged 
        
        # Store results
        Beta.out[[j]][m, ] <- as.vector(par_pdraws$beta)
        Sig.out[[j]][m]    <- par_pdraws$sig2
        Tau.out[[j]][m]    <- par_pdraws$tau
        Phi.out[[j]][m]    <- par_pdraws$phi
        
        # step 2: update imputed values (sample from preidctive distribution)
        pdraw_zj_imp <- rnorm(nrow(Z_mis), 
                        mean = (Z_mis %*% Beta.out[[j]][m, ]), 
                        sd = sqrt(Sig.out[[j]][m]))
      
        Imp.out[[j]][m, ] <- pdraw_zj_imp
        
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
}

