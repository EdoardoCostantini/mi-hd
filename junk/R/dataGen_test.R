### Title:    imputeHD-comp generating simple data to try out methods
### Author:   Edoardo Costantini
### Created:  2019-NOV-13
### Modified: 2019-NOV-13
### Notes:    generate data exclusively to get familiar with the methods
###           no simulation set up

# Functions for data generation
  #f = function(x1,x2,x3,x4,x5) {10*sin(pi*x1*x2) + 20*(x3-.5)**2 + 10*x4 + 5*x5}
  f = function(x1,x2,x3,x4) {1 + 10*x1 + 20*x2 + 30*x3 + 40*x4}
  pr = function(n,x) {rbinom(n, size = 1, prob = exp(x-2)/(1+exp(x-2))) }
  replace_na <- function(x, r) {replace(x, r==1, NA)}
  
  missDataGen <- function(n, p, corre = "n"){
  ## Description
    # The dependent variable depends on 5 variables (see f to understand how)
    # Missing values are forces on y and those 5 predictors based on the values
    # of a 6th variable. Hence, for this dataset, p >= 6
    # Input: 
    # - n: how many observations
    # - p: how many variables (>=)
    # Output: two datasets, one complete and one with missing at random datasets 
  ## Trial inputs
    # n = 100 # sample size
    # p = 6 # number of predictors (relevant + junk)
    # corre = c("n") #, "y") which type of correlation strucutre do you want n = none
    
    # Predictors
    Sigma = diag(p) # no correlation between predictors
    
    if (corre == "y"){ # random correlations
      p_coefs <- runif(p*p, -1, 1)
      s <- abs(diag(p)-1) * p_coefs
      s[lower.tri(s)] = t(s)[lower.tri(s)]
      diag(s) <- 1
      Sigma <- s
    } 
    X <- as.data.frame(mvtnorm::rmvnorm(n, 
                                        mean = rep(0, p), 
                                        sigma=Sigma ) )
    # Noise for DV
    Z <- rnorm(n, 0, 1)
    
    # Continuous DV
    y = f(X[,1],X[,2],X[,3],X[,4]) + Z
    
    # Dichotomous Responses
    y_dicho <- factor(pr(n, (5*X[,2]+3*X[,5])), labels = c("A", "B"))
    
    # Multinomial Responses (nominal)
    a <- c(0, 0, 0)
    b <- c(2, 4, 6)
    denom <- 1 + exp(a[1] + X[,3] * b[1]) +
                 exp(a[2] + X[,3] * b[2]) + 
                 exp(a[3] + X[,3] * b[3])
    vProb = cbind( 1/denom,             # baseline category
                   exp(a[1] + X[,3] * b[1])/denom, 
                   exp(a[2] + X[,3] * b[2])/denom, 
                   exp(a[3] + X[,3] * b[3])/denom )
    DV_location = t(apply(vProb, 1, rmultinom, n = 1, size = 1))
    DV <- apply(DV_location, 1, function(x) which(x==1))
    y_categ <- factor(DV, labels = c("A", "B", "C", "D"))
    
    # Multinomial Responses (ordinal, with multtinomial odered logit)
    u <- rlogis(nrow(X), location = 0, scale = 1)
    y_order <- X[,2] + u
    threshold <- c(-Inf,-1,0,1, Inf)
    y_order <- cut(y_order, threshold, ordered_result = TRUE)
    levels(y_order) <- c(1,2,3,4)
    
    # Missingness
    # Following a simple MAR mechanisms where the missingness on the
    # first 4 predictors, and the tree dependent variables depends on
    # the values of the 5th and 6th variable
    R <- as.data.frame(split(pr(n*9, X[, 5]+X[, 6]), 1:(4+5) )) # missingness for 3 Ys and first 5 Xs
    yI    <- replace(y, R[,1]==1, NA)
    y_dichoI <- replace(y_dicho, R[,2]==1, NA)
    y_categI <- replace(y_categ, R[,3]==1, NA)
    y_orderI <- replace(y_order, R[,4]==1, NA)
    Xinc <- X # initialize matrix for incomplete dataset
    for (p in 1:ncol(X[,2:6])) {
      Xinc[,p] <- replace(X[,p], R[,p]==1, NA)
    }
    return(list(complete = cbind(y_dicho, y_categ, y_order, y, X), 
                incomplete = cbind(y_dichoI, y_categI, y_orderI, yI, Xinc),
                complete_cont = cbind(y, X),
                incomplete_cont = cbind( yI, Xinc))
           )
  }
