# Project:   imputeHD-comp
# Objective: data generation functions
# Author:    Edoardo Costantini
# Created:   2020-06-23
# Modified:  2023-04-19
# Notes: 

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

  # If collinearity factor is not NA
  if (!is.null(cond$collinearity)) {
    if (!is.na(cond$collinearity)) {
      # MAR predictors within-block-1 correlation
      Sigma[c(4, 5), c(4, 5)] <- cond$collinearity
      # MAR predictors within-block-2 correlation
      Sigma[c(9, 10), c(9, 10)] <- cond$collinearity
      # Within-block-3 correlation
      Sigma[-c(parms$blck1, parms$blck2), -c(parms$blck1, parms$blck2)] <- cond$collinearity
    }
  }

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

# Experiment 2 - Latent Structure -----------------------------------------

simData_lv <- function(parms, cond){
  # cond <-  conds[1,]
  # parms$n <- 1e4
  n_it_tot <- parms$n_it * cond$lv
  
  # Structural parameters of the measurement model
  # Latent Variables Covariance matrix 
# (1) Block structure/Autoregressive correlation structure
  Phi <- diag(cond$lv)
  # Block 1: highly correlated variables
  Phi[parms$blck1, ] <- parms$blck1_r
  # Block 2: not so highly correlated variables
  Phi[parms$blck2, ] <- parms$blck2_r
  # Block 3: uncorrelated variables
  block3_p <- cond$lv - length(c(parms$blck1, parms$blck2))
  Phi[-c(parms$blck1, parms$blck2), ] <- .01
  # Fix diagonal
  diag(Phi) <- 1
  # Make symmetric
  Phi[upper.tri(Phi)] <- t(Phi)[upper.tri(Phi)]
  # Make it covariance instead of correlation matrix (depending of lv_var)
  Phi <- Phi * sqrt(parms$lv_var) * sqrt(parms$lv_var)
  
  # Factor loadings (random factor)
# (2) Sample from unif
  if(cond$fl == "none"){
    lambda <- rep(0, n_it_tot)
  } else {
    if(cond$fl == "high"){
      lambda <- runif(n_it_tot, .9, .97)
    } else {
      lambda <- runif(n_it_tot, .5, .6)
    }
  }
  
  # Observed Items Covariance matrix
# (3) uncorrelated observation errors
  # For decisions reagrding the paramter values look into your 
  # PhD_diary notes
  Theta <- diag(n_it_tot)
  for (i in 1:length(lambda)) {
    # Theta[i, i] <- 1 - lambda[i]^2 * 1
    Theta[i, i] <- parms$item_var - lambda[i]^2 * Phi[1,1]
  }
  
# (4) Items Factor Complexity = 1 (see Bollen1989 p234
#     aka simple measurement structure)
  Lambda <- matrix(nrow = n_it_tot, ncol = cond$lv)
  start <- 1
  for (j in 1:cond$lv) {
    end <- (start + parms$n_it) - 1
    vec <- rep(0, n_it_tot)
    vec[start:end] <- lambda[start:end]
    Lambda[, j] <- vec
    start <- end + 1
  }
  
  # Sample Scores
# (5) scores on latent variable and items errors centered around 0
  scs_lv    <- rmvnorm(parms$n, rep(0, cond$lv), Phi)
    colnames(scs_lv) <- paste0("lv", 1:ncol(scs_lv))
  scs_delta <- rmvnorm(parms$n, rep(0, n_it_tot), Theta)
    colnames(scs_delta) <- paste0("d_z", 1:ncol(scs_delta))
    
  # Compute Observed Scores
# (6) items
  x <- matrix(nrow = parms$n, ncol = n_it_tot)
  for(i in 1:parms$n){
    # x[i, ] <- t(Lambda %*% scs_lv[i, ] + scs_delta[i, ])
    x[i, ] <- t(parms$item_mean + Lambda %*% scs_lv[i, ] + scs_delta[i, ])
  }
  colnames(x) <- paste0("z", seq(1:n_it_tot))

  # Function output
  return( list(dat    = as.data.frame(x),
               Phi    = Phi,
               Theta  = Theta,
               Lambda = Lambda,
               scores_lv = scs_lv) )
}

imposeMiss_lv <- function(dat_in, parms, cond){
  ## Description
  # Given a fully observed dataset and a param object containing 
  # info on the (latent) imposition mechanism, it returns a version of
  # the original data with imposed missingness on all the items 
  # indicated as target int parms$z_m_id
  ## Example Inputs
  # cond <- conds[1,]
  # dat_in   <- simData_lv(parms, cond)
  
  # Body
  # Define non-response vector
  dat_out <- dat_in$dat
  for (i in 1:parms$zm_n) {
    nR <- simMissingness(pm    = cond$pm,
                         data  = dat_in$dat,
                         preds = parms$rm_x,
                         type  = parms$missType,
                         beta  = parms$auxWts)
    
    # Fill in NAs
    dat_out[nR, parms$z_m_id[i]] <- NA
  }
  
  # Result
  return( dat_out )
}

# Experiment 4 ------------------------------------------------------------

imposeMiss_evs <- function(dat_in, parms, cond){
  ## Description
  # Given a fully observed dataset and a param object containing the regression
  # coefficients of the model to impose missingness, it returns a version of
  # the original data with imposed missingness on all the variables indicated
  # as target int parms$z_m_id
  ## Example Inputs
  # cond <- conds[1,]
  # dat_in   <- data_source[sample(1:nrow(data_source),
  #                                cond$n,
  #                                replace = TRUE), ]
  
  # Body
  # Define non-response vector
  dat_out <- dat_in
  rm_x    <- dat_in[, parms$rm_x]
  
  # Recode percieved threat from immigrants 1 = low, 10 = high
  # rm_x[,1] <- match(rm_x[,1], max(rm_x[,1]):min(rm_x[,1]) )
  
  # Compute stuff
  for (i in 1:parms$zm_n) {
    nR <- simMissingness(pm    = runif(1, parms$pm[1], parms$pm[2]),
                         data  = rm_x,
                         preds = parms$rm_x,
                         type  = parms$missType,
                         beta  = parms$auxWts)
    
    # Fill in NAs
    dat_out[nR, parms$z_m_id[i]] <- NA
  }
  
  # Result
  return( dat_out )
}

# Function to fix factors in the EVS population generation procedure
fix.factor <- function(v, dt_h, dt_imp){
  if(is.factor(dt_imp[[v]])){
    original <- data.frame(value = dt_imp[[v]])
    key <- data.frame(value = val_labels(dt_h[[v]]),
                      label = names(val_labels(dt_h[[v]])),
                      row.names = NULL)
    v.out <- factor(plyr::join(original, key)[, 2])
  } else {
    v.out <- dt_imp[[v]]
  }
  return(v.out)
}
