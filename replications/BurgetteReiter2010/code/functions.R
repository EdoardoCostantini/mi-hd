### Title:    Replication Burgette Reiter 2010 - Functions
### Author:   Edoardo Costantini
### Created:  2020-JAN-10
### Modified: 2020-JAN-10
### Notes:    Create simulated data for sim study p.1072

gen_yhat_BR1073 <- function(x = X_test, b = pool_CART_est){
  outcome <- b[1] + # intercept
    b[2]*x[,1] + b[3]*x[,2] + b[4]*x[,3] + b[5]*x[,8] + b[6]*x[,9] + 
    b[7]*x[,3]**2 + b[8]*x[,1]*x[,2] + b[9]*x[,8]*x[,9]
  return(outcome)
}

desc_mis <- function(missData){
  # Given a dataset with missing values, it returns the vairablewise number
  # of missing values and the total percentage of observed cases (to check
  # the set up by Burgette Reiter 2010)
  
  output <- c(
    colMeans(is.na(missData)),
    mean(rowSums(is.na(missData)) == 0)
  )
  names(output)[length(output)] <- "tot"
  
  return(output)
  
}

# Extract results ---------------------------------------------------------

rMSE <- function (y_pred, y_obs) {
  # output: root mean squared (prediction) error
  MSE <- sqrt( mean((y_obs - y_pred)^2) )
  return(MSE)
}

rmse_est <- function(x) {
  # Applyes rmse to all columns of a dataset
  # x <- output_List[[2]][[1]]
  apply(x, 2, rMSE, y_obs = x[, ncol(x)])
}

bias_est <- function(x, x_true) {
  # returns the bias of an esitmate x
  x - x_true
} 

avg_bias <- function(out, meth_indx = 1, parms){
  # Given an output list coming from a parLapply of singleRun(), it returns the
  # average bias across all repetitions for a given imputation method (index
  # numerically)
    
  meth_colnames <- colnames(out[[1]]$pool_est)[meth_indx]
  B_true <- parms$b_true
  
  bias_out <- matrix(rep(NA, parms$b*parms$dt_rep), ncol = parms$dt_rep)
  for (dt in 1:length(out)) {
    bias_out[, dt] <- bias_est(out[[dt]]$pool_est[, meth_colnames], B_true)
  }

  return(round(rowMeans(bias_out), 3))
  
}

avg_ci_cov <- function(out, col_indx, parms){
  # Given an output list coming from a parLapply of singleRun(), it returns the
  # proportion of times the ci, pooled after an imputation done according to
  # a specific method (col_indx selects the column corresponding to the method
  # in the output), of each parameter contains the true parameter (specified in
  # params)
  
  B_true <- parms$b_true
  
  ci_out <- 0
  for (dt in 1:length(out)) {
    confi <- out[[dt]]$pool_conf[, col_indx]
    ci_out <- ci_out + as.numeric(c(B_true > confi[,1] & B_true < confi[,2]))
  }
  
  bwise_cov <- ci_out/length(out)
  tot_cov <- sum(ci_out)/(length(B_true)*length(out))
  
  return(c(bwise_cov, tot_cov))
  
}

avg_rmse <- function(out, col_indx, parms){
  # Given an output list coming from a parLapply of singleRun(), it returns the
  # mean value of the rmse obtained with the pooled coef estiamtes resulting
  # from an impuation method defined by the col_indx
  
  output <- rep(NA, parms$dt_rep)
  
  for (dt in 1:length(out)) {
    y_hat <- out[[dt]]$y_hat[, col_indx]
    y_obs <- out[[dt]]$y_test
    rMSE(y_hat, y_obs)
    output[dt] <- rMSE(y_hat, y_obs)
  }
  
  return(round(mean(output), 3))
}

avg_miss <- function(out, parms){
  # Given an output list coming from a parLapply of singleRun(), it returns the
  # mean proportion of missingness per variable and overall
  
  output <- matrix(rep(NA, parms$dt_rep*(parms$p+2)), nrow = parms$dt_rep)
  
  for (dt in 1:length(out)) {
    output[dt, ] <- out[[dt]]$miss_descrps
  }
  
  output <- colMeans(output)
  names(output) <- names(out[[1]]$miss_descrps)
  
  return(round(output, 3))
}

# CART methods specifics --------------------------------------------------

bbootstrap <- function(x) { # Bayesian Bootstrap
  # Input: a variable of any type (x)
  #   examples:
  #     @x <- rownames(airquality)
  #     @x <- rbinom(30, 1, .5) #another example
  # Output: a bayesian bootstrap sample of the size of x
  # Used in: CART_impute
  # Notes: based on Rubin 1998
  size <- length(x)
  u <- sort(c(runif(length(x)-1, 0, 1))) # n-1 uniform draws
  g <- numeric(0)
  for (i in 1:(length(x))) {
    if(length(u[i-1]) == 0) u_prev <- 0 else u_prev <- u[i-1]
    g[i] <- u[i]-(u_prev)
    if(i == length(x)) {
      u[i] <- 1
      g[i] <- 1-(u[i-1])
    }
    #print(cbind(u[i], u_prev, g[i]) ) # check that it works
  }
  #sum(g)
  bbsample <- sample(x, 
                     size = size, 
                     replace = TRUE, 
                     prob = g)
  return(bbsample)
}


# Parallel stuff ----------------------------------------------------------

## Broadcast the library function of a list of packages:
applyLib <- function(pkgList)
  lapply(pkgList, library, character.only = TRUE, logical = TRUE)
