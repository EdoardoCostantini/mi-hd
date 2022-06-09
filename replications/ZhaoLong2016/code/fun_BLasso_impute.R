### Title:    imputeHD-comp impute w/ Regularized Frequentiest Regressions DURR
### Author:   Anonymized for peer review
### Created:  2020-05-04
### Notes:    Function can be used only for univariate missingness imputation

impute_BLAS_hans <- function(Xy, Xy_mis, chains = 5, iter_bl = 1e3, burn_bl = 1e2){
  
  # Prep data ---------------------------------------------------------------
  # chains = parms$chains
  # iter_bl = parms$iter_bl
  # burn_bl = parms$burn_bl
  imputed_datasets <- vector("list", chains)
    names(imputed_datasets) <- seq(1, chains)
  
  j <- which(colSums(is.na(Xy_mis)) != 0)
  O <- 1 * (!is.na(Xy_mis))

  zj_obs <- as.vector(Xy_mis[O[, j] == 1, j])
  zj_mis <- as.vector(Xy_mis[O[, j] == 0, j])
  Z_obs <- as.matrix(Xy_mis[O[, j] == 1, -j])
  Z_mis <- as.matrix(Xy_mis[O[, j] == 0, -j])

  beta0 <- rep(0, ncol(Z_obs))

  # Sample parameters from posterior ----------------------------------------
  model <- blasso::blasso.vs(zj_obs, Z_obs,
                             iters = iter_bl,
                             burn = burn_bl,
                             beta = beta0,
                             beta.prior = "scaled",
                             sig2 = 1, sig2prior = c(a = .1, b = .1),
                             tau = 1, tauprior = c(r = .01, s = .01),
                             phi = .5, phiprior = c(h = 1, g = 1),
                             fixsig = FALSE, fixtau = FALSE, fixphi = FALSE)

  sample_indx <- sample(nrow(model$beta), chains)
  beta_dot <- model$beta[sample_indx, ]
  b_dot <- cbind(model$mu[sample_indx], beta_dot)

  sigma_dot_2 <- model$sig2[sample_indx]

  # predict
  zj_imp <- apply(b_dot, 1, function(b) {rnorm(nrow(Z_mis),
                                               mean = (Z_mis %*% b),
                                               sd = sqrt(sigma_dot_2))} )
  
  for(i in 1:chains){
    imputed_datasets[[i]] <- Xy_mis
    imputed_datasets[[i]][O[, j] == 0, j] <- zj_imp[, i]
  }
  
  return(imputed_datasets)
  
}
