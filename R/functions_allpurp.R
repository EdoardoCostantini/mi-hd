### Title:    imputeHD-comp all purpuse functions
### Author:   Edoardo Costantini
### Created:  2019-DEC-2
### Modified: 2019-DEC-2

# Functions ---------------------------------------------------------------

detect_family <- function(x){
# input: a dv for a glmnet::glmnet model 
#   examples:
#       @x = mtcars$mpg # gaussian
#       @x = iris$Species # multinomial
# output: the family that glmnet should use based on the DV type
# used in: DURR
  if(!is.null(levels(x))){
    family <- ifelse(length(levels(x))>2, "multinomial", "binomial")
  } else {
    family <- "gaussian" # limited to nomrally distributed var for now
  }
  return(family)
}

unigaus.lf <- function(theta,y,X){
# input: a set of model parameters theta, a dv y, and a set of predictors X
#   examples:
#     @theta <- c(1,1,1) # intercept, b1, sigma
#     @X <- cbind(1,runif(100))
#     @y <- X %*% theta[1:2] + rnorm(100, 0, theta[3])
# output: negative log likelihood of the data given some theta values
# used in: IURR for finding MLE of model parameters
  n <- nrow(X)
  k <- ncol(X)
  beta <- theta[1:k]
  sigma2 <- theta[k+1]
  e <- y - X%*%beta
  logl<- -.5*n*log(2*pi)-.5*n*log(sigma2)-((t(e)%*%e)/(2*sigma2))
  return(-logl)
}

arrange_dt_x_k <- function(dt, var_name){
# Input: a complete dataset, initialized or augmented, (dt) and the name the under-imputation variable
#   examples:
#     @dt       <- mtcars
#     @var_name <- names(mtcars)[3]
# Output: rearranged dt with the under imputation variable in the first spot (ready to be passed in a 
#         tree growing function)
# Used in: CART_impute
  indx_k <- names(dt) %in% var_name         # indexing obj: which variable is under imputation
  indx_no_k <- !names(dt) %in% var_name     # indexing obj: all vairables EXCEPT the one under imputation 
  Xno_k <- dt[, indx_no_k]                  # Use augmented Xno_k 
  x_k <- dt[, indx_k]                       # Select variable under imputation
  output <- data.frame(x_k = x_k, Xno_k)    # Combine the two for growing the tree
  # (doubt: should I grow the tree only on the observed values?)
  return(output)
}


tree_type <- function(x){ # Tree type definition
# Input: a measured variable that needs to be imputed (x)
#   examples:
#     @x <- mtcars$mpg            # continuous example
#     @x <- as.factor(mtcars$am)  # dichotomous example
# Output: "class" / "anova" the value of the method argument in the rpart tree growing function
#         (defines classification or a regression tree, respectively)
# Used in: CART_impute
  if(is.factor(x)) {type <- "class"} else {type <- "anova"}
  return(type)
}

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