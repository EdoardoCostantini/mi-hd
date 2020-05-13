### Title:    imputeHD-comp impute w/ Regularized Frequentiest Regressions DURR
### Author:   Edoardo Costantini
### Created:  2020-05-04
### Notes:    Function can be used only for univariate missingness imputation

impute_BLAS_hans <- function(Xy, Xy_mis, chains=5, iters=5, iter_bl = 1, burn_bl = 1e3, parms){
  
  # Prep data ---------------------------------------------------------------
  # chains = 2
  # iters = 2
  # iter_bl = 1
  # burn_bl = 1e3
  
  Z  <- Xy_mis
  p  <- ncol(Z) # number of variables [indexed with j]
  
  r  <- colSums(!is.na(Z)) # variable-wise count of observed values
  vmis_ind <- r != nrow(Z)
  p_imp <- sum(vmis_ind)
  
  # To store imputed values and check convergence
  imp_blasso_val <- vector("list", parms$chains)
    names(imp_blasso_val) <- seq(1:parms$chains)
  # empty storing objects for imputate values
  imps <- vector("list", p_imp)
    names(imps) <- names(which(r != nrow(Z)))
    for (v in 1:p_imp) {
      imps[[v]] <- matrix(rep(NA, (iters*(nrow(Z)-r[v]))), 
                          ncol = iters)
    }
  
  # To store multiply imputed datasets
  imp_blasso_dat <- vector("list", parms$chains)
    names(imp_blasso_dat) <- seq(1:parms$chains)
  
  O <- !is.na(Z)
  
  # Time performance
  start.time <- Sys.time()
  
  # For one chain of imputatuions
  for (cc in 1:chains) {
    
    Zm <- init_dt_i(Z, missing_type(Z)) # initialize data for each chain
    
    for (i in 1:iters) {
      print(paste0("BLASSO - Chain: ", cc, "/", chains, "; Iter: ", i, "/", iters))
      for (j in 1:p_imp) {
        zj_obs <- as.vector(Zm[O[, j] == TRUE, j])
        zj_mis <- as.vector(Zm[O[, j] == FALSE, j])
        Z_obs <- as.matrix(Zm[O[, j] == TRUE, -j])
        Z_mis <- as.matrix(Zm[O[, j] == FALSE, -j])
        
        beta0 <- rep(0, ncol(Z_obs)) # starting vbalues for gibbs sampler
        # step 1: take posterior sample for theta_hat_m_j
        model <- blasso::blasso.vs(Y = zj_obs, X = Z_obs,
                                   iters = iter_bl, 
                                   burn = burn_bl,
                                   beta = beta0, 
                                   beta.prior = "scaled", 
                                   sig2 = 1, sig2prior = c(a = .1, b = .1),
                                   tau = 1, tauprior = c(r = .01, s = .01), 
                                   phi = .5, phiprior = c(h = 1, g = 1),
                                   fixsig = FALSE, fixtau = FALSE, fixphi = FALSE,
                                   noisy = FALSE) # silent progress
        theta_hat <- as.vector(model$beta)
        sigma_hat_2 <- model$sig2
        
        # step 2: sample from preidctive distribution
        zj_imp <- rnorm(nrow(Z_mis), 
                        mean = (Z_mis %*% theta_hat), 
                        sd = sqrt(sigma_hat_2))
        # Append imputation
        imps[[j]][, i] <- zj_imp
        Zm[!O[, j], j] <- zj_imp
      }
    }
    
    imp_blasso_dat[[cc]] <- Zm
    imp_blasso_val[[cc]] <- imps
    
  }
  
  end.time <- Sys.time()
  
  return(list(dats = imp_blasso_dat,
              imps = imp_blasso_val,
              time = difftime(end.time, 
                              start.time, 
                              units = "mins"))
  )
}
