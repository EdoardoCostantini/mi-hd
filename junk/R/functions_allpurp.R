### Title:    imputeHD-comp all purpuse functions
### Author:   Edoardo Costantini
### Created:  2019-DEC-2
### Modified: 2020-JAN-9

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

unigaus.lf.dnorm<-function(theta,y,X){
# input: a set of model parameters theta, a dv y, and a set of predictors X
#   examples:
#     @theta <- c(1,1,1) # intercept, b1, sigma
#     @X <- cbind(1,runif(100))
#     @y <- X %*% theta[1:2] + rnorm(100, 0, theta[3])
# output: negative log likelihood of the data given some theta values
# used in: IURR for finding MLE of model parameters
  n<-nrow(X)
  k<-ncol(X)
  beta <- theta[1:k]
  sigma2 <- theta[k+1]
  e <- y - X%*%beta
  R = suppressWarnings(dnorm(e, 0, sqrt(sigma2), log=TRUE)) # less error prone
  return(-sum(R))
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

missing_type <- function(Z){
# Input: a dataset with missing values
#   examples:
#     @Z <- mice::boys
# Output: a list containing the names of the variables to be imputed in different formats
# Used in: MICERandomForest, and many more
# Notes: integer and numeric variables are considered (inaccurately) both as continuous
  l <- ncol(Z) - sum(colSums(is.na(Z)) == 0) # number of variables needing imputation
  l_names <- names(which(colSums(is.na(Z)) != 0)) # select the names of these k variables
  
  # Define names of variables w/ missing values (general and by measurment scale)
  vartypes <- rbind(lapply(lapply(Z[, names(Z) %in% l_names], class), paste, collapse = " "))
  contVars <- colnames(vartypes)[vartypes == "numeric" | vartypes == "integer"] # names of continuous variables for selection
  factVars <- colnames(vartypes)[vartypes == "factor"] # includes factors and ordered factors
  ordeVars <- colnames(vartypes)[vartypes == "ordered factor"] # includes factors and ordered factors
  
  output <- list(l=l, l_names=l_names, 
                 vartypes=vartypes, contVars=contVars, factVars=factVars, ordeVars=ordeVars)
  
  return(output)
}

init_dt_i <- function(Z0, missVarInfo){
  # Input: (1) a dataset with missing values; (2) and object produced by function missing_type
  #   examples:
  #     @Z <- mice::boys
  #     @missVarInfo <- missing_type(Z)
  # Output: a dataset with cotninuous variables imputed with mean value, and categorical with mode category
  # Used in: MICERandomForest
  # Notes: integer and numeric variables are considered (inaccurately) both as continuous
  
  # Make oredered factors as numeric
  Z0[, missVarInfo$ordeVars] <- as.numeric(Z0[, missVarInfo$ordeVars])
  
  # Impute sample means for continuous variables and odered factors
  s_means <- apply(Z0[, c(missVarInfo$contVars, missVarInfo$ordeVars)], 2, mean, na.rm = TRUE) # sample means
  for (j in 1:length(c(missVarInfo$contVars, missVarInfo$ordeVars))) {
    Z0 <- Z0 %>% mutate_at(vars(c(missVarInfo$contVars, missVarInfo$ordeVars)[j]),
                           ~replace(., is.na(.), s_means[j])
    )
  }
  
  # Impute most common level for unordered factors
  for (j in 1:length(missVarInfo$factVars)) {
    x <- addNA(Z0[, missVarInfo$factVars[j]])
    m_commo <- names(which.max(table(x)))
    levels(x) <- c(levels(Z0[, missVarInfo$factVars[j]]), 99)
    x[x == 99] <- m_commo
    x <- droplevels(x)
    Z0[, missVarInfo$factVars[j]] <- x
  }
  
  return(Z0)  # an initialized dataset
  # when dataset hass been itialized, m=0, so Z0 = {z0_1, z0_2, ... , z0_l, z_l+1, ... , z_p}
  # each z_j of this data will be the m-1 "previous iteration" version at the beginning of 
  # the variable loop (for j in 1:l) and the current iteration data at the end
}
