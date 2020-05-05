### Title:    imputeHD-comp impute w/ Regularized Frequentiest Regressions DURR
### Author:   Edoardo Costantini
### Created:  2020-02-26
### Modified: 2020-02-26
### Notes:    Function can be used only for univariate missingness imputation

impute_BLAS <- function(Xy_mis, chains = 5, iter_bl = 1e3, burn_bl = 1e2){
  
  # Prep data ---------------------------------------------------------------
  # chains = m
  # iter_bl = parms$iter_bl
  # burn_bl = parms$burn_bl
  imputed_datasets <- vector("list", chains)
  names(imputed_datasets) <- seq(1, chains)
  
  j <- which(colSums(is.na(Xy_mis)) != 0)
  O <- 1 * (!is.na(Xy_mis))
  
  y_obs <- as.vector(Xy_mis[O[, j] == 1, j])
  y_mis <- as.vector(Xy_mis[O[, j] == 0, j])
  X_obs <- as.matrix(Xy_mis[O[, j] == 1, -j])
  X_mis <- as.matrix(Xy_mis[O[, j] == 0, -j])
  X <- as.matrix(Xy_mis[, -j])
  
  # Sample parameters from posterior ----------------------------------------
  ############################### #
  ### Version monomvn::blasso ### #
  ############################### #
  model <- monomvn::blasso(X = X_obs, y = y_obs,
                      T = iter_bl,
                      # number of posterior samples to keep
                      thin = burn_bl,
                      RJ = TRUE,
                      rao.s2 = TRUE,
                      normalize = TRUE,
                      lambda2 = 0.1, # initial lasso penality
                      ab = c(.1, .1), # sigma**2 prior: IG(a,b)
                      rd = c(.01, .01), # lamba prior: Gamma(r,s)
                      mprior = c(1, 1), # rho prior: Beta(g,h)
                      verb = 1 #verbosity lvl
  )

  sample_indx <- sample(burn_bl : nrow(model$beta), chains)
  beta_dot <- model$beta[sample_indx, ]
  b_dot <- cbind(model$mu[sample_indx], beta_dot)

  sigma_dot_2 <- model$s2[sample_indx]
  # ############################### #
  
  # blOut <- bl(data       = as.data.frame(cbind(y_obs, X_obs)),
  #             y          = "y_obs",
  #             X          = setdiff(colnames(Xy_mis), c("y_obs")),
  #             iterations   = c(100, 100, 100),
  #             # Here thre numbers aprox, tuning sampling
  #             sampleSizes  = list(rep(100, 2),  # burn-in, retain for approximation 
  #                                 rep(1e3, 2), # burn-in, retain tuning
  #                                 rep(1e3, 2)) # burn-in, retain sampling bigger or same size
  #             # Approximates gets in the right neighbour, tuning improves, and sampling
  #             # 100 for burn in
  #             # Look into Gelman book for the thre approach
  # )
  # 
  # modelbeta <- getParams(blOut, 1)$beta
  # sample_indx <- sample(burn_bl : nrow(modelbeta), chains)
  # beta_dot <- modelbeta[sample_indx, ]
  # b_dot <- beta_dot
  # sigma_dot_2 <- getParams(blOut, 1)$sigma[sample_indx]
  
  y_dat <- apply(b_dot, 1, function(x) {rnorm(nrow(X_mis), 
                                              mean = (cbind(1, X) %*% x), 
                                              sd = sqrt(sigma_dot_2))} )
  for(i in 1:chains){
    imputed_datasets[[i]] <- Xy_mis
    imputed_datasets[[i]][O[, j] == 0, j] <- y_dat[, i]
  }
  
  return(imputed_datasets)
  
}
