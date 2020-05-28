### Title:    Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-05-19

# data generation ---------------------------------------------------------

genData <- function(cnd, parms){
  # For internals
  # cnd <- conds[1,]
  
  p_Zm <- length(parms$z_m_id) # number of variables w/ missing values
  p_Zf <- cnd["p"] - p_Zm      # number of fully observed variables
  
  # Gen Zf(ully observed)
  Zf <- rmvnorm(n     = parms$n, 
                mean  = rep(0, p_Zf), 
                sigma = AR1(p_Zf, rho = cnd["rho"]) )

  # Append empty Zm(issing values on p_Zm variables)
  Zm <- matrix(rep(NA, p_Zm*parms$n),
               ncol = p_Zm)
  Z <- cbind(Zm, Zf)
  
  # Fill Zm with true values
  AS_indx <- which(names(parms$S_all) == paste0("q", cnd["q"])) # active set (AS) indx
  Zs <- Z[, parms$S_all[[AS_indx]]] # Active set
  a <- rep(1, ncol(Zs)) * parms$stnr[AS_indx]
  
  Zm <- sapply(1:p_Zm, 
               function(x) {
                 rnorm(parms$n, 1 + Zs %*% a, sqrt(parms$z_m_var))
               }
  )
  
  # Replace in dataset Z the new true values
  Z <- cbind(Zm, Zf)
    colnames(Z) <- paste0("z", 1:ncol(Z))
  
  # Generate y
  b0 <- parms$b
  b <- rep(parms$b, parms$k)
  y <- rnorm(parms$n, 
             mean = b0 + Z[, 1:parms$k] %*% b,
             sd = sqrt(parms$y_var))
  
  return( as.data.frame(cbind(Z, y)) )
}

# Example use
# set.seed(1234)
# Xy <- genData(conds[1, ], parms)
# lm(y ~ z1 + z2 + z3 + z4 + z5, data = Xy) # true y model

imposeMiss <- function(Xy, parms){
  # Given a fully observed dataset and param object containing the regression
  # coefficients of the model to impose missingness, it returns a version of
  # the original data with imposed missingness on z1, and its missingness
  # idicator (R)
  
  n <- nrow(Xy)
  coefs <- parms$b_miss_model
  Xy_miss <- as.matrix(Xy)
  nR <- sapply(1:length(parms$z_m_id), function(d){
    logit_miss <- cbind(1, Xy_miss[, parms$detlamod[[d]]]) %*% coefs
    prb_miss <- exp(logit_miss)/(1+exp(logit_miss))
    nR <- rbinom(n, 1, prb_miss) == 1
    return(nR)
  }
  )
  
  for (i in 1:length(parms$z_m_id)) {
    Xy_miss[nR[, i], parms$z_m_id[i]] <- NA
  }
  
  return(list(Xy_miss = as.data.frame(Xy_miss),
              nR = nR) )
}

# generic functions -------------------------------------------------------

update_report <- function(method_name = "Method #", rep_count = 1, parms){
  # Updates the report .txt file with the count and name of accomplished task
  cat(paste0(Sys.time(), " | Reptetition ", rep_count, ": ", method_name, " done",
             "\n"),
      file = paste0(parms$outDir, parms$report_file_name),
      append = TRUE)  
}

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

missing_type <- function(Z){
  # Input: a dataset with missing values
  #   examples:
  #     @Z <- mice::boys
  # Output: a list containing the names of the variables to be imputed in different formats
  # Notes: - integer and numeric variables are considered (inaccurately) both as continuous;
  #        - dataframes are subsetted with [] to preserve structure
  
  l <- ncol(Z) - sum(colSums(is.na(Z)) == 0) # number of variables needing imputation
  l_names <- names(which(colSums(is.na(Z)) != 0)) # select the names of these k variables

  # By variable
  vartypes <- rbind(lapply(lapply(Z[colnames(Z) %in% l_names], 
                                  class), 
                           paste, collapse = " ")) # paste fixes ord factors
  
  # By measurement scale
  contVars <- colnames(vartypes)[vartypes == "numeric" | vartypes == "integer"]
    # names of continuous variables for selection
  factVars <- colnames(vartypes)[vartypes == "factor"] 
    # includes factors and ordered factors
  ordeVars <- colnames(vartypes)[vartypes == "ordered factor"] 
    # includes factors and ordered factors
  
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
  if( (length(missVarInfo$ordeVars))!=0 ){
    Z0[, missVarInfo$ordeVars] <- as.numeric(Z0[, missVarInfo$ordeVars])
  }

  # Impute sample means for continuous variables and odered factors
  if( (length(missVarInfo$contVars))!=0 ){
    s_means <- apply(Z0[c(missVarInfo$contVars, missVarInfo$ordeVars)], 
                     2, 
                     mean, 
                     na.rm = TRUE) # sample means
    
    for (j in 1:length(c(missVarInfo$contVars, missVarInfo$ordeVars))) {
      Z0 <- Z0 %>% mutate_at(vars(c(missVarInfo$contVars, missVarInfo$ordeVars)[j]),
                             ~replace(., is.na(.), s_means[j])
      )
    }
  }
  
  # Impute most common level for unordered factors
  if( (length(missVarInfo$factVars))!=0 ){
    for (j in 1:length(missVarInfo$factVars)) {
      
      x <- Z0[, missVarInfo$factVars[j]]
      m_commo <- names(which.max(table(x)))
      x[is.na(x)] <- m_commo
      Z0[, missVarInfo$factVars[j]] <- x
      
    }
  }
  
  return(Z0)  # an initialized dataset
  # when dataset has been itialized, m=0, so Z0 = {z0_1, z0_2, ... , z0_l, z_l+1, ... , z_p}
  # each z_j of this data will be the m-1 "previous iteration" version at the beginning of 
  # the variable loop (for j in 1:l) and the current iteration data at the end
}

# process_4IURR <- function(Z, Zm, j_th, parms){
#   ## Description:
#   # Given a datset with missing values, an imputed version from a previous
#   # iteration of imputation, and the index of the variable under imputation
#   # it returns a dataset with observed values on the variable under imputation
#   # and corresponding values on the other X variables either observed or 
#   # previously imputed.
#   ## For internals:
#   # data(iris2) # example dataset from PcAux package
#   # Z <- iris2[,-c(1,7)] # original data w/ missing values
#   # Zm <- init_dt_i(Z, missing_type(Z)) # result of previous iteration
#   # j_th <- 1
#   
#   ## Step 0. Prepare data for j-th variable imputation
#   ry <- !is.na(Z[, j_th]) # indexing obserrved values on y
#   Wm_j  <- Zm[,-j_th]     # predictors for imp model from inizialized dataset (m-1) [ALL CASES]  
#   zm_j  <- Zm[,j_th]      # outcome for imp model from inizialized dataset (m-1)    [ALL CASES]
#   Wm_mj <- Wm_j[!ry, ] # predictor rows for cases with missing z_j               [OBSERVED and IMPUTED]
#   zm_mj <- zm_j[!ry]   # current (m) draw of z_j for cases with missing z_j      [IMPUTED/MISSING CASES]
#   
#   # Select cases for imputation model
#   z_j_obs  <- Zm[, j_th][!is.na(Z[, j_th])]   # observed components of j-th variable [OBSERVED VALUES]
#   Wm_j_obs <- Zm[, -j_th][!is.na(Z[, j_th]),] # current (m) components of the predictors cooresponding 
#                                         # to observed cases on j-th variable [OBSERVED and IMPUTED VALUES]
#   
#   Zy <- (list(Wm_j_obs = as.matrix(Wm_j_obs),
#               z_j_obs = as.vector(z_j_obs),
#               Wm_mj = as.matrix(Wm_mj),
#               zm_mj = as.vector(zm_mj))) 
#   
#   return(Zy)
# }


# Estimation --------------------------------------------------------------

rr_est_lasso <- function(X, y, parms){
  ## Description:
  # Given any dv y (e.g. variable to be imputed), and its corresponding X 
  # values (observed, or imputed), fits a lasso regression
  ## For internals:
    # data("Boston", package = "MASS")
    # X <- model.matrix(medv~., Boston)
    # y <- Boston$medv
  ## Internals from simualtion
  # X = Wm_j_obs
  # y = z_j_obs
  # parms = parms
  
  cv_lasso <- cv.glmnet(X, y,
                        family = "gaussian",
                        nfolds = 10,
                        alpha = 1)
  
  regu.mod <- glmnet(X, y, 
                     family = "gaussian",
                     alpha = 1, 
                     lambda = cv_lasso$lambda.min)
  
  return(regu.mod)
}

rr_est_elanet <- function(X, y, parms){
  ## Description:
  # Given any dv y (e.g. variable to be imputed), and its corresponding X 
  # values (observed, or imputed), fits a lasso regression
  ## For internals:
  #   data("Boston", package = "MASS")
  #   X <- model.matrix(medv~., Boston)
  #   y <- Boston$medv
  # Set training control
  train_control <- trainControl(method = "cv",
                                number = 10,
                                # repeats = 5,
                                selectionFunction = "best",
                                verboseIter = FALSE)
  
  # CV
  el_cv <- train(y ~., 
                 data = cbind(y, X), 
                 method = "glmnet",
                 family = "gaussian", 
                 type.multinomial = "grouped",
                 trControl = train_control,
                 preProcess = c("center", "scale"),
                 tuneLength = 25
  )
  
  # Train model
  regu.mod <- glmnet(X, y,
                     family = "gaussian",
                     type.multinomial = "grouped", # if glmfam is multinomial, otherwise does not apply
                     alpha = el_cv$bestTune[1],
                     lambda = el_cv$bestTune[2])
  
  return(regu.mod)
}


# Imputation --------------------------------------------------------------

imp_gaus_DURR <- function(model, X_tr, y_tr, X_te, parms){
  ## Description ##
  
  # Given a fitted imputation model it returns the imputations for
  # the missing values on a normally distributed imputation dv
  
  ## For internals ##
  
  # data("Boston", package = "MASS")
  # train_ind <- Boston$medv %>%
  #   createDataPartition(p = 0.8, list = FALSE)
  # train  <- Boston[train_ind, ]
  # test <- Boston[-train_ind, ]
  # 
  # X_tr <- model.matrix(medv~., train)[,-1]
  # X_te <- model.matrix(medv~., test)[,-1]
  # y_tr <- train$medv
  #   cv <- cv.glmnet(X_tr, y_tr, alpha = 1) # cross validate lambda value
  # model <- glmnet(X_tr, y_tr, alpha = 1, lambda = cv$lambda.min)
  
  ## Body ##
  s2hat   <- mean((predict(model, X_tr) - y_tr)**2) 
    # according to paper this is the estimate
  X_te    <- model.matrix(rep(1, nrow(X_te)) ~., data.frame(X_te))[, -1]
    # create a prediction matrix (for possible factor cooding needed),
    # -1 no need for intercept
  yhat_te <- rnorm(n =    nrow(X_te),
                   mean = predict(model, X_te), 
                   sd =   sqrt(s2hat))
  return(yhat_te)
}

imp_dich_DURR <- function(model, X_tr, y_tr, X_te, parms){
  ## Description ##
  
  # Given a fitted imputation model it returns the imputations for
  # the missing values on dichotomous
  
  ## For internals ##
  
  # dt_dich <- data.frame(x = matrix(rnorm(100 * 20), 100, 20),
  #                       g2 = as.factor(sample(0:1, 100, replace = TRUE)))
  # train_ind <- dt_dich$g2 %>%
  #   createDataPartition(p = 0.8, list = FALSE)
  # train  <- dt_dich[train_ind, ]
  # test <- dt_dich[-train_ind, ]
  # 
  # X_tr <- model.matrix(g2~., train)[,-1]
  # X_te <- model.matrix(g2~., test)[,-1]
  # y_tr <- train$g2
  #   cv <- cv.glmnet(X_tr, y_tr, family = "binomial", alpha = 1)
  # model <- glmnet(X_tr, y_tr,
  #                 family = "binomial", alpha = 1, lambda = cv$lambda.min)
  
  ## Body ##
  
  X_te <- model.matrix(rep(1, nrow(X_te)) ~ X_te)[, -1]
  py1  <- predict(model, X_te, type = "response") # probability y = 1
  y_prd <- rbinom(nrow(X_te), 1, py1)             # random draw
  
  # Fixing outcome levels
  y_prd <- factor(y_prd, levels(y_tr)) # return to original levels
  
  return(y_prd)
}

imp_multi_DURR <- function(model, X_tr, y_tr, X_te, parms){
  ## Description ##
  
  # Given a fitted imputation model it returns the imputations for
  # the missing values on a categorical variable w/ multiple classes
  
  ## For internals ##
  
  # dt_multi <- data.frame(x = matrix(rnorm(100 * 20), 100, 20),
  #                       g4 = as.factor(sample(1:4, 100, replace = TRUE)))
  # train_ind <- dt_multi$g4 %>%
  #   createDataPartition(p = 0.8, list = FALSE)
  # train  <- dt_multi[train_ind, ]
  # test <- dt_multi[-train_ind, ]
  # 
  # X_tr <- model.matrix(g4~., train)[,-1]
  # X_te <- model.matrix(g4~., test)[,-1]
  # y_tr <- train$g4
  #   cv <- cv.glmnet(X_tr, y_tr, family = "multinomial", alpha = 1)
  # model <- glmnet(X_tr, y_tr,
  #                 family = "multinomial", alpha = 1, lambda = cv$lambda.min)

  ## Body ##
  X_te <- model.matrix(rep(1, nrow(X_te)) ~ X_te)[, -1]
  py   <- predict(model, X_te, type = "response")
  DV_location = t(apply(py, 1, rmultinom, n = 1, size = 1)) # 1 for the category in which is most
  # likely that an observation is
  colnames(DV_location) <- levels(y_tr)
  z.m_j_mis <- apply(DV_location, 1, function(x) names( which(x==1) ))
  
  return(z.m_j_mis)
}

imp_gaus_IURR <- function(model, X_tr, y_tr, X_te, y_te, parms){
  ## Description ##
  
  # Given a regularized model it returns the imputations for
  # the missing values on a normally distributed imputation dv
  # according to IURR method
  
  ## For internals ##

  # data("Boston", package = "MASS")
  # train_ind <- Boston$medv %>%
  #   caret::createDataPartition(p = 0.8, list = FALSE)
  # train  <- Boston[train_ind, ]
  # test <- Boston[-train_ind, ]
  # 
  # X_tr <- model.matrix(medv~., train)[,-1]
  # X_te <- model.matrix(medv~., test)[,-1]
  # y_tr <- train$medv
  # y_te <- test$medv
  #   cv <- cv.glmnet(X_tr, y_tr, alpha = 1) # cross validate lambda value
  # model <- glmnet(X_tr, y_tr, alpha = 1, lambda = cv$lambda.min)
  
  # model = regu.mod
  # X_tr = Wm_j_obs
  # y_tr = z_j_obs
  # X_te = Wm_mj 
  # y_te = zm_mj
  # parms = parms
  
  ## Body ##
  
  # Select predictors based on rr
    rr_coef <- as.matrix(coef(model)) # regularized regression coefs
    rr_coef_no0 <- row.names(rr_coef)[rr_coef != 0]
    AS <- rr_coef_no0[-1] # predictors active set
  
  # MLE estimate of model parameters
    
  # 1. define starting values
    if(identical(rr_coef_no0, "(Intercept)")){
      lm_fit <- lm(y_tr ~ 1)
      X_mle <- model.matrix(y_tr ~ 1)
    } else {
      X_mle <- model.matrix(y_tr ~ X_tr[, AS])
      lm_fit <- lm(y_tr ~ X_tr[, AS])
    }
    startV <- c(coef(lm_fit), sigma(lm_fit))
    
  # 2. optimize loss function
    MLE_fit <- optim(startV, 
                     .lm_loss,
                     method="BFGS",
                     hessian=T,
                     y = y_tr, X = X_mle)
    
  # 3. obtain estimates
    theta <- MLE_fit$par
    OI <- solve(MLE_fit$hessian) # parameters cov maatrix
    
  # Sample parameters for posterior predictive distribution
    pdraws_par <- MASS::mvrnorm(1, 
                                mu = MLE_fit$par, 
                                Sigma = OI)
    
  # Sample posterior predictive distribution
    if(identical(rr_coef_no0, "(Intercept)")){
      y_imp <- rnorm(n = nrow(X_te),
                     mean = pdraws_par[1],
                     sd = pdraws_par[2])
    } else {
      X_ppd <- model.matrix( ~ X_te[, AS]) # X for posterior pred dist
      b_ppd <- pdraws_par[-length(pdraws_par)] # betas for posterior pred dist
      sigma_ppd <- tail(pdraws_par, 1) # sigma for posterior pred dist
      y_imp <- rnorm(n = nrow(X_te),
                     mean = X_ppd %*% b_ppd,
                     sd = sigma_ppd)
    }
  return(y_imp)
}

.lm_loss <-function(theta, y, X){
  n <- nrow(X)
  k <- ncol(X)
  beta <- theta[1:k]
  sigma <- theta[k+1]
  if(sigma < 0) {
    dev <- 10000000
  } else {
    ll <- dnorm(y, mean = X %*% beta, sd = sigma, log = TRUE)
    dev <- -2 * sum(ll)
  }
  # Return 
  return(dev)
}

.lambda <- function(x) (x + 1) / (x + 3)

.miDf <- function(m, b, t, dfCom) {
  fmi   <- .fmi(m, b, t)
  df0   <- (m - 1) * (1 / fmi^2)
  dfObs <- .lambda(dfCom) * dfCom * (1 - fmi)
  
  df0 / (1 + (df0 / dfObs))
}

.fmi <- function(m, b, t){
  # proportion of variation attributable to the missing data
  # aka fmi (not adjusted for the finite number of imps)
  fmi <- (1 + 1/m) * b/t
  return(fmi)
}

fit_models <- function(multi_dt, mod){
  # Given a list of complete datasets it fits a model described
  # in mod 
  models <- lapply(X = multi_dt,
                   FUN = function(x) lm(mod, data = x))
  return(models)
}

get_pool_EST <- function(fits){
  ## Description
  # Given a list of imputed datasets under the same imputation model
  # it returns the pooled estiamtes of the regression coefs
  ## For internals
  # fake_data <- nhanes
  # colnames(fake_data) <- c("y", "z1", "z2", "z3")
  # imp <- mice(fake_data, print = FALSE, m = 10, seed = 24415)
  # imp_dats <- mice::complete(imp, "all")
  # multi_dt = imp_dats
  summa_models <- lapply(X = fits,
                         FUN = function(x) summary(x))
  coefs <- t(sapply(X = summa_models,
                    FUN = function(x) coef(x)[, "Estimate"]))
  Q_bar <- colMeans(coefs)
  return(Q_bar)
}

get_pool_CI <- function(fits){
  ## Description
  # Given a list of imputed datasets under the same imputation model
  # it returns the pooled CIs of the regression coefs
  ## For internals
  # fake_data <- nhanes
  # colnames(fake_data) <- c("y", "z1", "z2", "z3")
  # imp <- mice(fake_data, print = FALSE, m = 10, seed = 24415)
  # imp_dats <- mice::complete(imp, "all")
  # multi_dt = imp_dats
  m <- length(fits)
  
  summa_models <- lapply(X = fits,
                         FUN = function(x) summary(x))
  ## Coef estimates ##
  coefs <- t(sapply(X = summa_models,
                    FUN = function(x) coef(x)[, "Estimate"]))
  Q_bar <- colMeans(coefs)
  
  ## Variances
  all_vcov <- lapply(X = summa_models,
                     FUN = function(x) vcov(x))
  
  U_bar <- diag(Reduce('+', all_vcov) / m)
  
  B <- diag(1 / (m-1) * (t(coefs) - Q_bar) %*% t(t(coefs) - Q_bar))
  
  T <- U_bar + B + B/m
  
  ## Degrees of freedom
  nu <- .miDf(length(fits), b = B, t = T, summa_models[[1]]$df[2])
 
  ## CI computation
  t_nu <- qt(1 - (1-parms$alphaCI)/2, 
             df = nu)
  
  CI <- data.frame(lwr = Q_bar - t_nu * sqrt(T), 
                   upr = Q_bar + t_nu * sqrt(T))
  return(CI = CI)
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

# Results -----------------------------------------------------------------

bias_est <- function(x, x_true) {
  # returns the bias of an esitmate x
  x - x_true
} 

check_cover <- function(x){ # 1 is the same value for all parameters
  # given a 2 x (p+1) matrix containing lower and uper CI bounds
  # it checks the shared parameter value 1 is included or not in
  # the interval
  return(x[, 1] < 1 & x[, 2] > 1)
}

extract_results <- function(cond_name, output, dt_rep){
  # Example input
  # cond_name <- names(out[[1]])[1]
  # output <- out
  # dt_rep = out[[1]]$cond_200_4$parms$dt_rep
  
  # Bias
  store_sum <- vector("list", dt_rep)
  
  for (i in 1:dt_rep) {
    store_sum[[i]] <- output[[i]][[cond_name]]$cond_bias
  }
  
  bias_out <- round(Reduce("+", store_sum)/dt_rep, 3)
  # bias_b1 <- as.data.frame(t(bias_out))[2] # only interested in b1
  bias <- as.data.frame(t(bias_out))#[1:4]
  
  # Average Coverage
  store_sum <- vector("list", dt_rep)
  
  for (i in 1:dt_rep) {
    store_sum[[i]] <- output[[i]][[cond_name]]$cond_CIco
  }
  
  CI_out <- Reduce("+", store_sum)/dt_rep
    rownames(CI_out) <- rownames(bias_out)
  # CI_b1 <- as.data.frame(t(CI_out))[2]
  CI <- as.data.frame(t(CI_out))#[1:4]
  
  # resu <- cbind(bias_b1, CI_b1)
  # colnames(resu) <- c("bias", "ci")
  resu <- list(bias = bias, 
               CI = round(CI, 3))
  return(resu)
}
