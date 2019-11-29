### Title:    imputeHD-comp generating simple data to try out methods
### Author:   Edoardo Costantini
### Created:  2019-NOV-13
### Modified: 2019-NOV-13
### Notes:    generate data exclusively to get familiar with the methods
###           no simulation set up

# Functions for data generation
  f = function(x1,x2,x3,x4,x5) {10*sin(pi*x1*x2) + 20*(x3-.5)**2 + 10*x4 + 5*x5}
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
    # Output: two datasets, one compelte and one with missing at random datasets 
  ## Trial inputs
    # n = 20 # sample size
    # p = 6 # number of predictors (relevant + junk)
    # corre = c("n") #, "y") which type of correlatio n strucutre do you want
    
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
    
    # Dependent
    # Model based on an example (Firedman) reported by Xu et al 2016 for the BART paper
    y = f(X[,1],X[,2],X[,3],X[,4],X[,5]) + Z # continuous
    
    # Dichotomous Responses
    y_dicho <- factor(pr(n, f(X[,1],X[,2],X[,3],X[,4],X[,5])), labels = c("A", "B"))
    
    # Multinomila Responses (nominal)
    a <- c(4, 2, 1)
    b <- c(2, 4, 6)
    denom <- 1 + exp(a[1] + X[,1] * b[1]) +
                 exp(a[2] + X[,1] * b[2]) + 
                 exp(a[3] + X[,1] * b[3])
    vProb = cbind( 1/denom,             # baseline category
                   exp(a[1] + X[,1] * b[1])/denom, 
                   exp(a[2] + X[,1] * b[2])/denom, 
                   exp(a[3] + X[,1] * b[3])/denom )
    DV_location = t(apply(vProb, 1, rmultinom, n = 1, size = 1))
    DV <- apply(DV_location, 1, function(x) which(x==1))
    y_categ <- factor(DV, labels = c("A", "B", "C", "D"))
    
    # Multinomila Responses (ordinal, with multtinomial odered logit)
    u <- rlogis(nrow(X), location = 0, scale = 1)
    y <- X[,2] + u
    threshold <- c(-Inf,-1,0,1, Inf)
    y_order <- cut(y, threshold, ordered_result = TRUE)
    levels(y_order) <- c(1,2,3,4)
    
    # Missingness
    # Following a simple MAR mechanisms where the missingness on the
    # first 5 predictors, and the tree dependent variables depends on
    # the values of the 6th variable
    R <- as.data.frame(split(pr(n*9, X[, 6]), 1:(4+5) )) # missingness for 3 Ys and first 5 Xs
    yI    <- replace(y, R[,1]==1, NA)
    y_dichoI <- replace(y_dicho, R[,2]==1, NA)
    y_categI <- replace(y_categ, R[,3]==1, NA)
    y_orderI <- replace(y_order, R[,4]==1, NA)
    Xinc <- X # initialize matrix for incomplete dataset
    for (p in 1:ncol(X[,2:6])) {
      Xinc[,p] <- replace(X[,p], R[,p]==1, NA)
    }
    return(list(complete = cbind(y, X), 
                incomplete = cbind(y_dichoI, y_categI, y_orderI, yI, Xinc)))
  }

# # Do it out of function (to save the data)
# # dimensionality
#   n = 100 # sample size
#   p = 200 # number of predictors (relevant + junk)
#   
# 
# # Predictors
#   set.seed(20191113)
#   X <- as.data.frame(mvtnorm::rmvnorm(n, 
#                                       mean = rep(0, p), 
#                                       sigma=diag(rep(1, p))) )
#   round(apply(X, 2, var),1)
#   round(cor(X), 2) # just at random you get some correlation between certail predictors
#   
# # Noise
#   Z <- rnorm(n, 0, 1)
# 
# # Dependent
#   # Model based on an example (Firedman) reported by Xu et al 2016 for the BART paper
#   y = f(X[,1],X[,2],X[,3],X[,4],X[,5]) + Z
# 
# # Missingness
#   # Following a simple MAR mechanisms where the missingness on the
#   # first 5 predictors depends on the values of the 6th variable
#   R <- as.data.frame(split(rep( pr(n, X[, 6]), 6), 1:(1+ncol(X[,1:5])) )) # missingness for Y and first 5 Xs
#   yinc <- replace(y, R[,1]==1, NA)
#   Xinc <- X # initialize matrix for incomplete dataset
#   for (p in 1:ncol(X[,2:6])) {
#     Xinc[,p] <- replace(X[,p], R[,p]==1, NA)
#   }
# 
# # Missing data type
#   sum(rowSums(is.na(cbind(yinc, Xinc[, 1:5]))) != 0)/n # percetnage of cases with missing values 
#   mice::md.pattern(cbind(yinc, Xinc[, 1:5]))           # missing data patterns
#   mice::md.pairs(cbind(yinc, Xinc[, 1:5]))
# 
# # Save data
#   out <- list(complete = cbind(y, X), 
#               incomplete = cbind(yinc, Xinc))
#   saveRDS(out,"../data/data_tryout.rds")
  