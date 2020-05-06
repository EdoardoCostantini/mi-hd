### Title:    test IURR approach
### Author:   Edoardo Costantini
### Created:  2020-03-05
### Modified: 2020-03-05
### Notes:    test what might go wrong with the IURR impute function

source("./init.R")

cond <- conds[1,]
m <- parms$chains
reg_type = "lasso"

for (i in 1:500) {
  Xy <- genData(cond, parms)
  Xy_mis <- imposeMiss(Xy, parms)$Xy_miss
  
  # imp_IURR_lasso <- impute_IURR(Xy_mis = Xy_mis, 
  #                               cond = cond, 
  #                               chains = m,
  #                               reg_type = "lasso")
  
  j <- which(colSums(is.na(Xy_mis)) != 0)
  O <- 1 * (!is.na(Xy_mis))
  
  z_j_obs  <- y_obs <- as.vector(Xy_mis[O[, j] == 1, j])
  zm_mj    <- y_mis <- as.vector(Xy_mis[O[, j] == 0, j])
  Wm_j_obs <- X_obs <- as.matrix(Xy_mis[O[, j] == 1, -j])
  Wm_mj    <- X_mis <- as.matrix(Xy_mis[O[, j] == 0, -j])
  X <- as.matrix(Xy_mis[, -j])
  
  
  # Fit regularized regression ----------------------------------------------
  # Not elegant way of accepting only regularized fits if the selected 
  # predictors are less than the number of observations
  p_check <- ncol(Wm_j_obs)
  n_check <- nrow(Wm_j_obs)
  
  while (p_check > n_check) {
    
  # Lasso
  if(reg_type == "lasso"){
    regu.mod <- rr_est_lasso(X = Wm_j_obs, y = z_j_obs, parms = parms)
  }
  # Elastic net
  if(reg_type == "el"){
    regu.mod <- rr_est_elanet(X = Wstar_mj, y = zstar_j, parms = parms)
  }
    
  model <- regu.mod
  X_tr = Wm_j_obs
  y_tr = z_j_obs
  X_te = Wm_mj
  y_te = zm_mj
  
  lasso.coef <- as.matrix(coef(model))
  S_hat <- row.names(lasso.coef)[lasso.coef != 0][-1] # estimated active set
  
  # Estimate Imputation model on observed z_j
  lm_fit <- lm(y_tr ~ X_tr[, S_hat])
  
  theta_MLE <- coef(lm_fit)
  Sigma_MLE <- vcov(lm_fit)
  
  # Check status
  solve(vcov(lm_fit))
  p_check <- sum(as.vector(regu.mod$beta) != 0)
  n_check <- regu.mod$nobs # nrow(Wm_j_obs)
  
  }
  
  # Impute ------------------------------------------------------------------
  # zm_j <- imp_gaus_IURR(model = regu.mod,
  #                       X_tr = Wm_j_obs, y_tr = z_j_obs,
  #                       X_te = Wm_mj, y_te = zm_mj, 
  #                       parms = parms)
  
  # Sample parameters for prediction/imputation
  theta_m_j <- MASS::mvrnorm(1, 
                             theta_MLE, 
                             Sigma_MLE + diag(ncol(Sigma_MLE)) * parms$k_IURR )
  sigma_m_j <- sigma(lm_fit)
  # This is not correct. I'm using the OLS estiamtes with no variation in
  # sigma. Later fix this!
  
  # Get imputations
  y_imp <- rnorm(n = nrow(X_te),
                 mean = (model.matrix( ~ X_te[, S_hat]) %*% theta_m_j),
                 sd = sigma_m_j)
  print(i) 
}
