### Title:    Replication Burgette Reiter 2010 - Functions
### Author:   Edoardo Costantini
### Created:  2020-JAN-10
### Modified: 2020-JAN-10
### Notes:    Create simulated data for sim study p.1072

gen_X_BR1073 <- function(n = 1e3){
  # input: a dv for a glmnet::glmnet model 
  #   examples:
  #    n <- 1e3 # sample size
  #    b0 <- 0 # intercept
  #    B <- c(0, .5,.5,.5,.5,.5,1,1) # regression coefficient values
  #    eps <- rnorm(n)
  # output: complete dataset based on simulation set up by Burgette Reiter 2010
  
  # gen X
  # Covariance
  C1 <- matrix(rep(.5, 4*4), ncol = 4)
    diag(C1) <- 1
    C1[lower.tri(C1)] <- NA
  C2 <- matrix(rep(.3, 10*10), ncol = 10)
    diag(C2) <- 1
    C2[lower.tri(C2)] <- NA
    C2[1:nrow(C1), 1:ncol(C1)] <- C1
  
  tC2 <- t(C2)
  tC2[upper.tri(C2)] <- C2[upper.tri(C2)]
  
  Sigma_X <- diag(rep(1, 10)) %*% tC2 %*% diag(rep(1, 10)) 
    # assuming 1 as sd for every variable
  mu_X <- rep(0, 10)
  
  X <- rmvnorm(n, mu_X, Sigma_X)
  
  return(X)
}

gen_y_BR1073 <- function(x, 
                         b0 = 0, 
                         b = c(0, .5,.5,.5,.5,.5,1,1)){
  # input: a dv for a glmnet::glmnet model 
  #   examples:
  #    x <- mtcars[, 1:7]
  #    b0 <- 0 # intercept
  #    B <- runif(ncol(x)) # regression coefficient values
  # output: dependent variable for B&R2010 set up
  
  y <- b0 + 
    b[1]*x[,1] + b[2]*x[,2] + b[3]*x[,3] + b[4]*x[,8] + b[5]*x[,9] + 
    b[6]*x[,3]**2 + b[7]*x[,1]*x[,2] + b[8]*x[,8]*x[,9] + 
    rnorm(nrow(x))
  
  return(y)
}

imp_miss_BR1073 <- function(){
  
}

gen_yhat_BR1073 <- function(x = X_test, b = pool_CART_est){
  outcome <- b[1] + # intercept
    b[2]*x[,1] + b[3]*x[,2] + b[4]*x[,3] + b[5]*x[,8] + b[6]*x[,9] + 
    b[7]*x[,3]**2 + b[8]*x[,1]*x[,2] + b[9]*x[,8]*x[,9]
  return(outcome)
}

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

