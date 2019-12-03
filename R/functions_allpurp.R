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
