### Title:    imputeHD-comp impute w/ Regularized Frequentiest Regressions DURR
### Author:   Edoardo Costantini
### Created:  2020-05-04
### Notes:    Function can be used only for univariate missingness imputation

impute_BLAS_hans <- function(Xy, Xy_mis, m = 5, iters = 10, iter_bl = 1e3, burn_bl = 1e2){
  
  # Prep data ---------------------------------------------------------------
  # chains = m
  # iter_bl = 1
  # burn_bl = 1e3

  # TRIAL SCALE
  Xy <- genData(cond, parms)
  Xy$z1 <- Xy$z1+100
  Xy_mis <- imposeMiss(Xy, parms)$Xy_miss
  Xy[1:3]
  Xy_mis[1:3]
  
  Z  <- Xy_mis
  p  <- ncol(Z) # number of variables [indexed with j]
  
  r  <- colSums(!is.na(Z)) # variable-wise count of observed values
  vmis_ind <- r != nrow(Z)
  p_imp <- sum(vmis_ind)
  
  # To store imputed values and check convergence
  imputed_values <- vector("list", p_imp)
  names(imputed_values) <- names(which(r != nrow(Z)))
  for (v in 1:p_imp) {
    imputed_values[[v]] <- matrix(rep(NA, (iters*(nrow(Z)-r[v]))), 
                                  ncol = iters)
  }
  
  # To store multiply imputed datasets
  imputed_datasets <- vector("list", iters)
  names(imputed_datasets) <- seq(1, iters)
  
  O <- !is.na(Xy_mis)
  
  # For one chain of imputatuions
  Zm <- init_dt_i(Xy_mis, missing_type(Xy_mis)) # initialize data for each chain
  for (j in 1:p_imp) {
    j <- 1
    # Scale data for blasso before selecting variables and based on original data
    # so that you can revert to the orignal scale the final imputations
    zj <- Zm[, j]
    Z <- Zm[ , -j]
    
    round(colMeans(Z), 1)
    sapply(Z, var)
    
    round(colMeans(Zm), 1)
    sapply(Zm, var)
    
    Xy_mis_scale <- Xy_mis
    Xy_mis_scale[, j] <- zj
    Xy_mis_scale[, -j] <- Z
    
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
                               fixsig = FALSE, fixtau = FALSE, fixphi = FALSE)
    theta_hat <- as.vector(model$beta)
    sigma_hat_2 <- model$sig2
    
    # step 2: sample from preidctive distribution
    # scale Z_mis according to what blaso did to Z_obs
    sZ_mis <- sapply(1:ncol(Z_mis), function(x) {
      Z_mis[, x] <- (Z_mis[, x]-mean(Z_obs[, x]))/sd(Z_obs[, x])
    })
    # predict on this "blasso" scale
    zj_imp <- rnorm(nrow(Z_mis), 
                    mean = (sZ_mis %*% theta_hat), 
                    sd = sqrt(sigma_hat_2))
    
    # Rescale zj_imp to original data
    zj_imp <- zj_imp + mean(zj_obs)
      # when running blasso, zj_obs is centered around its mean. Hence, these
      # predicted values need to be "rescaled" summing that mean so that
      # they are again on the orignal data scale
  }
  
  # Get required number of samples
  sample_indx <- sample(burn_bl : nrow(model$beta), chains)
  beta_dot <- model$beta[sample_indx, ]
  b_dot <- cbind(model$mu[sample_indx], beta_dot)

  sigma_dot_2 <- model$sig2[sample_indx]
  
  zj_imp <- apply(b_dot, 1, function(b) {rnorm(nrow(Z_mis), 
                                              mean = (Z_mis %*% b), 
                                              sd = sqrt(sigma_dot_2))} )
  # Rescale zj_imp to original data
  zj_imp <- zj_imp + mean(zj_obs)
    # when running blasso, zj_obs is centered around its mean. Hence, these
    # predicted values need to be "rescaled" summing that mean so that
    # they are again on the orignal data scale
  
  for(i in 1:chains){
    imputed_datasets[[i]] <- Xy_mis
    imputed_datasets[[i]][O[, j] == 0, j] <- zj_imp[, i]
  }
  
  return(imputed_datasets)
  
}
