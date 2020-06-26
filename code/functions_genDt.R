### Title:    data geneartion functions
### Porject:  Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-06-23

# data generation ---------------------------------------------------------

# multivairate 

genDt_mvn <- function(cond, parms){
  # For internals
  # cond <- conds[1,]
  
  # Generate covaraince matrix
  # For paramters decisions, look back at
  Sigma <- diag(cond$p)
  
  # Block 1: highly correlated variables
  Sigma[parms$blck1, ] <- parms$blck1_r
  
  # Block 2: not so highly correlated variables
  Sigma[parms$blck2, ] <- parms$blck2_r
  
  # Block 3: uncorrelated variables
  block3_p <- cond$p - length(c(parms$blck1, parms$blck2))
  Sigma[-c(parms$blck1, parms$blck2), ] <- .01 #runif(block3_p*cond$p, 0.01, .1)

  # Fix diagonal
  diag(Sigma) <- 1
  
  # Make symmetric
  Sigma[upper.tri(Sigma)] <- t(Sigma)[upper.tri(Sigma)]
  
  # Gen Axuliariy Variables
  Z <- rmvnorm(n     = parms$n, 
               mean  = rep(0, cond$p), 
               sigma = Sigma )
  colnames(Z) <- paste0("z", 1:ncol(Z))

  # Scale according to your evs examination
  Z_sc <- Z * sqrt(parms$item_mean) + parms$item_var
  # values used based on study of democratic and morals items
  
  return( Z_sc )
}

# Use
# Xy <- genDt_mvn(conds[1, ], parms)
# cor(Xy)[1:15, 1:15]
# cov(Xy)[1:15, 1:15]

imposeMAR <- function(target_id, dt_in, parms, cond){
  # Right Tail Imposition of missingness
  # Inputs
  # target_id <- parms$z_m_id[4]
  # dt_in   <- genDt_mvn(conds[1, ], parms)
  # cond <- conds[1,]
  
  # Body
  dt_out <- as.matrix(dt_in)
  rm_id <- which(parms$z_m_id %in% target_id)
  
  # Perform
  parms$rm_b <- c(-3, -3, 2, -1)
  linPred <- dt_out[, parms$rm_x[rm_id, ]] %*% parms$rm_b
  missPropensity <- pnorm(linPred, mean(linPred), sd(linPred))
  nR <- missPropensity >= 1 - cond$pm
  dt_out[nR, target_id] <- NA
  
  # Result
  return( dt_out[, target_id] )
}

imposeMiss <- function(dt_in, parms, cond){
  # Given a fully observed dataset and param object containing the regression
  # coefficients of the model to impose missingness, it returns a version of
  # the original data with imposed missingness on z1, and its missingness
  # idicator (R)
  
  ## Inputs ##
  
  # dt_in   <- genDt_mvn(conds[1, ], parms)
  # cond <- conds
  
  ## Body ##
  # n <- nrow(Xy)
  dt_out <- as.matrix(dt_in)

  for (i in parms$z_m_id) {
    dt_out[, i] <- imposeMAR(target_id = i,
                             dt_in = dt_in,
                             parms = parms,
                             cond = cond)
  }
  
  return( as.data.frame(dt_out) )
}

# DengEtAl Style ----------------------------------------------------------
# 
# genData <- function(cnd, parms){
#   # For internals
#   # cnd <- conds[1,]
#   
#   p_Zm <- length(parms$z_m_id) # number of variables w/ missing values
#   p_Zf <- cnd["p"] - p_Zm      # number of fully observed variables
#   
#   # Gen Zf(ully observed)
#   Zf <- rmvnorm(n     = parms$n, 
#                 mean  = rep(0, p_Zf), 
#                 sigma = AR1(p_Zf, rho = cnd["rho"]) )
#   
#   # Append empty Zm(issing values on p_Zm variables)
#   Zm <- matrix(rep(NA, p_Zm*parms$n),
#                ncol = p_Zm)
#   Z <- cbind(Zm, Zf)
#   
#   # Fill Zm with true values
#   AS_indx <- which(names(parms$S_all) == paste0("q", cnd["q"])) # active set (AS) indx
#   Zs <- Z[, parms$S_all[[AS_indx]]] # Active set
#   a <- rep(1, ncol(Zs)) * parms$stnr[AS_indx]
#   
#   Zm <- sapply(1:p_Zm, 
#                function(x) {
#                  rnorm(parms$n, 1 + Zs %*% a, sqrt(parms$z_m_var))
#                }
#   )
#   
#   # Replace in dataset Z the new true values
#   Z <- cbind(Zm, Zf)
#   colnames(Z) <- paste0("z", 1:ncol(Z))
#   
#   # Generate y
#   b0 <- parms$b
#   b <- rep(parms$b, parms$k)
#   y <- rnorm(parms$n, 
#              mean = b0 + Z[, 1:parms$k] %*% b,
#              sd = sqrt(parms$y_var))
#   
#   return( as.data.frame(cbind(Z, y)) )
# }
# 
# # Example use
# # set.seed(1234)
# # Xy <- genData(conds[1, ], parms)
# # lm(y ~ z1 + z2 + z3 + z4 + z5, data = Xy) # true y model
# 
# imposeMiss <- function(Xy, parms){
#   # Given a fully observed dataset and param object containing the regression
#   # coefficients of the model to impose missingness, it returns a version of
#   # the original data with imposed missingness on z1, and its missingness
#   # idicator (R)
#   
#   n <- nrow(Xy)
#   coefs <- parms$b_miss_model
#   Xy_miss <- as.matrix(Xy)
#   nR <- sapply(1:length(parms$z_m_id), function(d){
#     logit_miss <- cbind(1, Xy_miss[, parms$detlamod[[d]]]) %*% coefs
#     prb_miss <- exp(logit_miss)/(1+exp(logit_miss))
#     nR <- rbinom(n, 1, prb_miss) == 1
#     return(nR)
#   }
#   )
#   
#   for (i in 1:length(parms$z_m_id)) {
#     Xy_miss[nR[, i], parms$z_m_id[i]] <- NA
#   }
#   
#   return(list(Xy_miss = as.data.frame(Xy_miss),
#               nR = nR) )
# }
# 
# genDataCat <- function(cnd, parms){
#   # For internals
#   # cnd <- conds[1,]
#   
#   p_Zm <- length(parms$z_m_id) # number of variables w/ missing values
#   p_Zf <- cnd["p"] - p_Zm      # number of fully observed variables
#   
#   # Gen Zf(ully observed)
#   Zf <- rmvnorm(n     = parms$n, 
#                 mean  = rep(0, p_Zf), 
#                 sigma = AR1(p_Zf, rho = cnd["rho"]) )
#   
#   # Append empty Zm(issing values on p_Zm variables)
#   Zm <- matrix(rep(NA, p_Zm*parms$n),
#                ncol = p_Zm)
#   Z <- cbind(Zm, Zf)
#   
#   # Fill Zm with true values
#   AS_indx <- which(names(parms$S_all) == paste0("q", cnd["q"])) # active set (AS) indx
#   Zs <- Z[, parms$S_all[[AS_indx]]] # Active set
#   a <- rep(1, ncol(Zs)) * parms$stnr[AS_indx]
#   
#   # Number of vairables of each type
#   nvar_count <- 2
#   nvar_dicho <- 1
#   nvar_polyt <- 0
#   
#   Zm_count <- sapply(1:nvar_count, 
#                      function(x) {
#                        rnorm(parms$n, 1 + Zs %*% a, sqrt(parms$z_m_var))
#                      }
#   )
#   
#   Zm_dicho <- data.frame(
#     lapply(1:nvar_dicho, # lapply keeps them as factors
#            function(x) {
#              logit <- 1 + Zs %*% a
#              prob <- exp(logit)/(1+exp(logit))
#              y <- factor(rbinom(parms$n, 1, prob))
#            }
#     )
#   )
#   
#   # Zm_polyt <- data.frame(
#   #   lapply(1:nvar_polyt, # lapply keeps them as factors
#   #          function(x) {genMultiCat(K = 3, X = Zs)$y})
#   # )
#   
#   Zm <- data.frame(Zm_count, Zm_dicho) #, Zm_polyt)
#   
#   # Replace in dataset Z the new true values
#   Z <- cbind(Zm, Zf)
#   colnames(Z) <- paste0("z", 1:ncol(Z)) # fix names
#   
#   # Generate y
#   X <- as.matrix(cbind(Z[, 1:2],
#                        model.matrix( ~ Z[,3])[, -1],
#                        # model.matrix( ~ Z[,3])[, -1],
#                        Z[, -c(1:3)]))
#   
#   b0 <- parms$b
#   b <- rep(parms$b, parms$k)
#   y <- rnorm(parms$n, 
#              mean = b0 + X[, 1:parms$k] %*% b,
#              sd = sqrt(parms$y_var))
#   
#   return( as.data.frame(cbind(Z, y)) )
# }
# 
# # Example use
# # set.seed(1234)
# # Xy <- genDataCat(conds[1, ], parms)
# # lm(y ~ z1 + z2 + z3 + z4 + z5, data = Xy) # true y model
# 
# imposeMissCat <- function(Xy, parms){
#   # Given a fully observed dataset and param object containing the regression
#   # coefficients of the model to impose missingness, it returns a version of
#   # the original data with imposed missingness on z1, and its missingness
#   # idicator (R)
#   
#   n <- nrow(Xy)
#   coefs <- parms$b_miss_model
#   nR <- sapply(1:length(parms$z_m_id), function(d){
#     logit_miss <- as.matrix(cbind(1, Xy[, parms$detlamod[[d]]])) %*% coefs
#     prb_miss <- exp(logit_miss)/(1+exp(logit_miss))
#     nR <- rbinom(n, 1, prb_miss) == 1
#     return(nR)
#   }
#   )
#   
#   for (i in 1:length(parms$z_m_id)) {
#     Xy[nR[, i], parms$z_m_id[i]] <- NA
#   }
#   
#   return(list(Xy_miss = Xy,
#               nR = nR) )
# }
# 
# genMultiCat <- function(K = 3, X) {
#   # Data gen from multinomial logit model
#   # Notes: Randomly generated parameters, in the future I will choose them
#   # Example inputs
#   # K = 3 # number of categories for outcome vairable
#   # X = matrix(rnorm(1e3), ncol = 4) # some dataset
#   
#   # Body
#   a <- runif(K-1, 0, 5) # intercept terms
#   B <- matrix(runif((K-1)*ncol(X), 0, 5),
#               ncol = ncol(X)) # logit reg coefs
#   
#   # Denominator
#   denom_pt1 <- sapply(1:(K-1), 
#                       function(k) exp(a[k] + X %*% B[k, ]) )
#   denom <- 1 + rowSums(denom_pt1)
#   
#   # Calculating the matrix of probabilities for each category of the DV
#   vProb = cbind( 1/denom,             # baseline category
#                  sapply(1:(K-1), 
#                         function(k) exp(a[k] + X %*% B[k, ])/denom))
#   # each row has the probabilities for that specific individual of falling
#   # in each category of Y (i.e. each row sums to 1)
#   # rowSums(vProb)
#   
#   # Sample category for each individual
#   location = t(apply(vProb, 1, rmultinom, n = 1, size = 1)) 
#   
#   # Store result
#   y <- apply(location, 1, function(x) which(x==1))
#   y <- factor(y, labels = sample(words, K)) 
#   
#   return(list(y = y,
#               true_par = cbind(a, B))
#   )
# }