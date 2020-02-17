### Title:    Replication Burgette Reiter 2010 - Simulation Subroutines
### Author:   Edoardo Costantini
### Created:  2020-JAN-10
### Modified: 2020-JAN-10
### Notes:    Routines for the simulation study


# Generate 1 dataset Xy
gen_Xy <- function(n = 1e3, b0 = 0, b = c(0, .5,.5,.5,.5,.5,1,1)){
  # input: a dv for a glmnet::glmnet model 
  #   examples:
  #    n <- 1e3 # sample size
  #    b0 <- 0 # intercept
  #    b <- c(0, .5,.5,.5,.5,.5,1,1) # regression coefficient values
  # output: fully observed dataset based on simulation set up by B&R 2010
  
  # Correlation matrix
  C1 <- matrix(rep(.5, 4*4), ncol = 4)
  diag(C1) <- 1
  C1[lower.tri(C1)] <- NA
  C2 <- matrix(rep(.3, 10*10), ncol = 10)
  diag(C2) <- 1
  C2[lower.tri(C2)] <- NA
  C2[1:nrow(C1), 1:ncol(C1)] <- C1
  tC2 <- t(C2)
  tC2[upper.tri(C2)] <- C2[upper.tri(C2)]
  # this is the correlation matrix
  
  Sigma_X <- diag(rep(1, 10)) %*% tC2 %*% diag(rep(1, 10)) 
  # assuming 1 as sd for every variable
  
  mu_X <- rep(0, ncol(Sigma_X))
  
  X <- as.data.frame(rmvnorm(n, mu_X, Sigma_X))
  
  y <- b0 + 
    b[1]*X[,1] + b[2]*X[,2] + b[3]*X[,3] + b[4]*X[,8] + b[5]*X[,9] + 
    b[6]*X[,3]**2 + b[7]*X[,1]*X[,2] + b[8]*X[,8]*X[,9] + 
    rnorm(nrow(X))
  
  return(cbind(X, y))
  
}

# Impose missingess on 1 Xy dataset
imposeMiss <- function(Xy){
  # input: dataset coming from gen_Xy
  #   examples:
  #    Xy <- gen_Xy() # see gen function in this script
  # output: dataset w/ missing values based on set up by B&R 2010
  
  preds <- c("V8", "V9")
  col_indx <- colnames(Xy)[!(colnames(Xy) %in% preds)]
  
  type <- sample(c("high", "low", "center", "tails"), 1)
  
  for (mis in 1:length(col_indx)) {
    Xy[simMissingness(pm = .17, data = Xy,
                      preds = preds,
                      type  = type, beta  = runif(2) ), col_indx[mis]] <- NA
    # simMissingness comes from Kyle's script and simulates a nonresponse vector
    # according to the proportion of missings desired on a variable, with certain
    # predictors (causes) coming from a given datasets
  }
  
  return(Xy)
  
}

singleRun <- function(rp = 10, chains = 10, iters = 1, parms){
  

  ## Set up seed ----------------------------------------------------------- ##
  # Setup the PRNG (needs rlecuyer package):
  
  .lec.SetPackageSeed(rep(parms$seed, 6))
  if(!rp %in% .lec.GetStreams())
    .lec.CreateStream(c(1 : parms$nStreams))
  .lec.CurrentStream(rp)
  
  ## Data ------------------------------------------------------------------ ##
  # Gen 1 dataset w/ missing values; 1 fully-obs data for out-of-sample rmse
  
  Xy <- gen_Xy()
  Xy_mis <- imposeMiss(Xy)
  miss_descrps <- desc_mis(Xy_mis)
  Xy_test <- gen_Xy(n = 500)
  X_test <- Xy_test[, -ncol(Xy_test)]
  y_test <- Xy_test[, ncol(Xy_test)]
  
  ## Imputation ------------------------------------------------------------ ##
  # Impute m times the data w/ missing values w/ different methods
  
  # MICE cart Burgette Reiter
  imp_CART_bb <- miceImpHDv::mice(Xy_mis, m = chains, maxit = iters,
                                  meth = 'cart.bb', minbucket = 5)

  # MICE default
  imp_norm <- mice::mice(Xy_mis, m = chains, maxit = iters, meth = 'norm')
  
  # MICE cart
  imp_CART <- mice(Xy_mis, m = chains, maxit = iters,
                                  meth = 'cart', minbucket = 5)
  
  ## Pooling --------------------------------------------------------------- ##
  # Analyse and pool estimates of the m imputations for each method
  
  # MICE cart Burgette Reiter
  fit <- with(imp_CART_bb, 
              expr = lm(y ~ V1 + V2 + V3 + V8 + V9 + I(V3^2) + V1:V2 + V8:V9))

  pool_CART_est <- mice::pool(fit)$pooled[,1]
  
  pool_CART_conf <- summary(mice::pool(fit), 
                            conf.int = TRUE)[,c("2.5 %", "97.5 %")]
  
  # MICE default
  fit_norm <- with(imp_norm, 
                   expr = lm(y ~ V1 + V2 + V3 + V8 + V9 + I(V3^2) + V1:V2 + V8:V9))
  pool_norm_est <- mice::pool(fit_norm)$pooled[,1]
  
  pool_norm_conf <- summary(mice::pool(fit_norm), 
                            conf.int = TRUE)[,c("2.5 %", "97.5 %")]
  
  # MICE cart Burgette Reiter
  fit <- with(imp_CART, 
              expr = lm(y ~ V1 + V2 + V3 + V8 + V9 + I(V3^2) + V1:V2 + V8:V9))
  
  pool_CARTd_est <- mice::pool(fit)$pooled[,1]
  
  pool_CARTd_conf <- summary(mice::pool(fit),
                             conf.int = TRUE)[,c("2.5 %", "97.5 %")]
  
  # collect estiamtes
  pool_est <- cbind(pool_CART_est, pool_norm_est, pool_CARTd_est)
  pool_conf <- cbind(pool_CART_conf, pool_norm_conf, pool_CARTd_conf)
  
  ## Prediction ------------------------------------------------------------ ##
  # Predict outcome variable with pooled estimates (for out-of-sample rmse)
  
  y_hat <- apply(pool_est, 2, gen_yhat_BR1073, x = X_test)
  
  ## Store output ---------------------------------------------------------- ##
  
  output <- list(pool_est = pool_est,
                 pool_conf = pool_conf,
                 y_hat = y_hat,
                 y_test = y_test,
                 miss_descrps = miss_descrps)
  
  return(output)
}

