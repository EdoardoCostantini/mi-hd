### Title:    data geneartion functions
### Porject:  Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-06-23

# data generation ---------------------------------------------------------

# Experiment 1: Multivariate Data -----------------------------------------
simData_exp1 <- function(cond, parms){
  # For internals
  # cond <- conds[1,]
  
  # 1. Generate covariance matrix -----------------------------------------
    Sigma <- diag(cond$p)
  # Block 1: highly correlated variables
    Sigma[parms$blck1, ] <- parms$blck1_r
  # Block 2: not so highly correlated variables
    Sigma[parms$blck2, ] <- parms$blck2_r
  # Block 3: uncorrelated variables
    block3_p <- cond$p - length(c(parms$blck1, parms$blck2))
    Sigma[-c(parms$blck1, parms$blck2), ] <- .01
  # Fix diagonal
    diag(Sigma) <- 1
  # Make symmetric
    Sigma[upper.tri(Sigma)] <- t(Sigma)[upper.tri(Sigma)]
  
  # 2. Gen n x p data -----------------------------------------------------
    # Sample
    Z <- rmvnorm(n     = parms$n, 
                 mean  = rep(0, cond$p), 
                 sigma = Sigma )
    # Give meaningful names
    colnames(Z) <- paste0("z", 1:ncol(Z))
  
  # 3. Scale according to your evs examination
    Z_sc <- Z * sqrt(parms$item_var) + parms$item_mean

  return( as.data.frame( Z_sc ) )
}

# Use
# Xy <- simData_exp1(conds[1, ], parms)
# cor(Xy)[1:15, 1:15]
# cov(Xy)[1:15, 1:15]

# Experiment 2: Interactions ----------------------------------------------

simData_exp2 <- function(parms, inter = TRUE){
  # For internals
  # cond <- conds[1,]
  # inter = TRUE
  
  # Generate covaraince matrix
  # For paramters decisions, look back at
  Sigma <- diag(parms$p - length(parms$blcky))
  
  # Block 1: highly correlated variables
  blck1_indx <- parms$blck1 - length(parms$blcky)
  Sigma[blck1_indx, ] <- parms$blck1_r
  
  # Block 2: not so highly correlated variables
  blck2_indx <- parms$blck2-length(parms$blcky)
  Sigma[blck2_indx, ] <- parms$blck2_r
  
  # Block 3: uncorrelated variables
  blck3_indx <- (1:ncol(Sigma))[-c(blck1_indx, blck2_indx)]
  Sigma[blck3_indx, ] <- .01
  
  # Fix diagonal
  diag(Sigma) <- 1
  
  # Make symmetric
  Sigma[upper.tri(Sigma)] <- t(Sigma)[upper.tri(Sigma)]
  
  # Gen Axuliariy Variables
  Z <- rmvnorm(n     = parms$n, 
               mean  = rep(0, ncol(Sigma)), 
               sigma = Sigma )
  colnames(Z) <- paste0("z", 1:ncol(Z))
  
  # Gen y variables
  # Main Effects
  blcky_indx <- c(blck1_indx[1], blck2_indx[1], blck3_indx[1])
  
  # Interaction Effects
  inter_indx <- c(blcky_indx[1], blcky_indx[3])
  int_term <- apply(Z[, inter_indx], 1, prod)
  
  # Linear model predictors
  if(inter == TRUE){
    betas <- c(.333, -.333, -.1, .666)
  } else {
    betas <- c(.333+.222, -.333-.222, -.1-.222, 0)
  }
  
  # Signal to noise ratio
  Z_pred <- cbind(Z[, blcky_indx], int_term)
  signal <- t(parms$beta) %*% cov(Z_pred) %*% parms$beta
  sY     <- (signal / parms$r2) - signal
  
  omega2  <-  as.vector( parms$stnr * betas %*% cov(Z_pred) %*% betas )
  eps     <-  rnorm(parms$n, mean = 0, sd = sqrt(sY))
  
  # Create y
  y <- cbind(Z[, blcky_indx], int_term) %*% betas + eps
  
  # Combine data
  yX <- data.frame(y = y, Z)
  
  # Scale according to your evs examination
  yX_sc <- yX * sqrt(parms$item_var) + parms$item_mean

  return(yX_sc)
}

# Xy <- genDt_exp2(parms, inter = TRUE)
# lm_out <- lm(y ~ -1 + z1 + z11 + z21 + z1:z21, data = Xy)
# summary(lm_out)
# lm_out_sd <- lm(y ~ -1 + z1 + z11 + z21 + z1:z21, data = as.data.frame(scale(Xy)))
# summary(lm_out_sd)
# 
# Xy <- genDt_exp2(parms, inter = FALSE)
# lm(y ~ -1 + z1 + z11 + z21, data = Xy)
# lm(y ~ -1 + z1 + z11 + z21, data = as.data.frame(scale(Xy)))


# Missingness -------------------------------------------------------------

imposeMiss <- function(dat_in, parms, cond){
  ## Description
  # Given a fully observed dataset and a param object containing the regression
  # coefficients of the model to impose missingness, it returns a version of
  # the original data with imposed missingness on all the variables indicated
  # as target int parms$z_m_id
  ## Example Inputs
  # cond <- conds[1,]
  # dat_in   <- simData_exp1(cond, parms)
  
  # Body
  # Define non-response vector
  dat_out <- dat_in
  for (i in 1:parms$zm_n) {
    nR <- simMissingness(pm    = cond$pm,
                         data  = dat_in,
                         preds = parms$rm_x[i, ],
                         type  = parms$missType,
                         beta  = parms$auxWts)
    
    # Fill in NAs
    dat_out[nR, parms$z_m_id[i]] <- NA
  }
  
  # Result
  return( dat_out )
}