### Title:    Replication Zhao Long 2016 - Functions
### Author:   Edoardo Costantini
### Created:  2020-02-20
### Modified: 2020-02-20

# data generation ---------------------------------------------------------

genData <- function(par_conds, parms){
  # par_conds <- conds[1,] # example for par_conds
  
  # Extract condition parameters
  p   <- par_conds[[1]]
  rho <- par_conds[[2]]
  q   <- par_conds[[3]]
  
  # Extract fixed parameters
  n     <- parms$n
  b     <- rep(parms$b, (4))
  y_var <- parms$y_var
  ActSet<- which(names(parms$S_all) == paste0("q", q)) # active set (AS) indx
  stnr  <- parms$stnr[ActSet] # seòecion of signal-to-noise for given AS
  S     <- parms$S_all[[ActSet]] # variables for AS (predictors)
  
  # Gen Z-1
  Z_m1 <- rmvnorm(n, rep(0, (p-1)), AR1((p-1), rho))
  
  # Gen z1
  a0 <- 1 # Based on specification of DengEtAl because Zhao was not clear
  a <- rep(1, length(S)) * stnr
  z1 <- rnorm(n, 
              mean = a0 * Z_m1[, S] %*% a, 
              sd = sqrt(1))
  # lm(z1 ~ -1 + Z_m1[, S])
  # Gen Z
  Z <- cbind(z1, Z_m1)
  colnames(Z) <- paste0("z", 1:ncol(Z))
  Z_desing <- cbind(rep(1, nrow(Z)), Z)
  
  # Gen Y
  y <- rnorm(n, Z_desing[, 1:4] %*% b, sqrt(y_var))
  
  return( as.data.frame(cbind(Z, y)) )
}

imposeMiss <- function(Xy, parms){
  # Given a fully observed dataset and param object containing the regression
  # coefficients of the model to impose missingness, it returns a version of
  # the original data with imposed missingness on z1, and its missingness
  # idicator (R)
  
  n <- nrow(Xy)
  coefs <- parms$b_miss_model
  var_mismech <- c("z2", "z3", "y")
  
  Xy <- as.matrix(Xy)
  logit_miss <- cbind(1, Xy[, var_mismech]) %*% coefs
  prb_miss <- exp(logit_miss)/(1+exp(logit_miss))
  
  nR <- rbinom(n, 1, prb_miss) == 1
  
  Xy[nR, "z1"] <- NA
  
  return(list(Xy_miss = as.data.frame(Xy),
              nR = nR) )
}

# generic functions -------------------------------------------------------

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
  # when dataset hass been itialized, m=0, so Z0 = {z0_1, z0_2, ... , z0_l, z_l+1, ... , z_p}
  # each z_j of this data will be the m-1 "previous iteration" version at the beginning of 
  # the variable loop (for j in 1:l) and the current iteration data at the end
}

process_4DURR <- function(Z, Zm, j_th, parms){
  ## Description:
  # Given a datset with missing values, an imputed version from a previous
  # iteration of imputation, and the index of the variable under imputation
  # it returns a dataset with observed values on the variable under imputation
  # and corresponding values on the other X variables either observed or 
  # previously imputed.
  ## For internals:
  # data(iris2) # example dataset from PcAux package
  # Z <- iris2[,-c(1,7)] # original data w/ missing values
  # Zm <- init_dt_i(Z, missing_type(Z)) # result of previous iteration
  # j_th <- 1

  ry <- !is.na(Z[, j_th]) # indexing observed values on y
  
  # For estimation of imputation model
  Wm_j  <- Zm[,-j_th]
  zm_j  <- Zm[, j_th]
  
  # For prediction based on imputation model
  Wm_mj <- Wm_j[!ry, ]  # Xs for j-th var w/ originally missing values
  zm_mj <- zm_j[!ry]    # imputation of j-th var at i-1 iteration
  
  # Generate bootstrap sample
  indx_boSample <- sample(1:nrow(Wm_j), nrow(Wm_j), replace = TRUE)
  W.star.m_j    <- Wm_j[indx_boSample, ]
  z.star_j      <- zm_j[indx_boSample]
  
  # but select only real observed values out of such sample
  z.star_j_obs     <- z.star_j[!is.na(Z[indx_boSample, j_th])]         # OBSERVED values on j-th
  W.star.m_j_obs   <- W.star.m_j[!is.na(Z[indx_boSample, j_th]), ] # Xs for j-th var originally observed values [OBSERVED and i-1 IMPUTATIONS]
  
  # Data for regularized models
  Zx <- model.matrix(z.star_j_obs~., data.frame(z.star_j_obs, W.star.m_j_obs))[,-1]
  zy <- z.star_j_obs
  
  Zy <- (list(Wstar_mj = Zx,
              zstar_j = zy,
              Wm_mj = Wm_mj,
              zm_mj = zm_mj)) 
  
  return(Zy)
}

process_4IURR <- function(Z, Zm, j_th, parms){
  ## Description:
  # Given a datset with missing values, an imputed version from a previous
  # iteration of imputation, and the index of the variable under imputation
  # it returns a dataset with observed values on the variable under imputation
  # and corresponding values on the other X variables either observed or 
  # previously imputed.
  ## For internals:
  # data(iris2) # example dataset from PcAux package
  # Z <- iris2[,-c(1,7)] # original data w/ missing values
  # Zm <- init_dt_i(Z, missing_type(Z)) # result of previous iteration
  # j_th <- 1
  
  ## Step 0. Prepare data for j-th variable imputation
  ry <- !is.na(Z[, j_th]) # indexing obserrved values on y
  Wm_j  <- Zm[,-j_th]     # predictors for imp model from inizialized dataset (m-1) [ALL CASES]  
  zm_j  <- Zm[,j_th]      # outcome for imp model from inizialized dataset (m-1)    [ALL CASES]
  Wm_mj <- Wm_j[!ry, ] # predictor rows for cases with missing z_j               [OBSERVED and IMPUTED]
  zm_mj <- zm_j[!ry]   # current (m) draw of z_j for cases with missing z_j      [IMPUTED/MISSING CASES]
  
  # Select cases for imputation model
  z_j_obs  <- Zm[, j_th][!is.na(Z[, j_th])]   # observed components of j-th variable [OBSERVED VALUES]
  Wm_j_obs <- Zm[, -j_th][!is.na(Z[, j_th]),] # current (m) components of the predictors cooresponding 
                                        # to observed cases on j-th variable [OBSERVED and IMPUTED VALUES]
  
  Zy <- (list(Wm_j_obs = as.matrix(Wm_j_obs),
              z_j_obs = as.vector(z_j_obs),
              Wm_mj = as.matrix(Wm_mj),
              zm_mj = as.vector(zm_mj))) 
  
  return(Zy)
}


# Estimation --------------------------------------------------------------



rr_est_lasso <- function(X, y, parms){
  ## Description:
  # Given any dv y (e.g. variable to be imputed), and its corresponding X 
  # values (observed, or imputed), fits a lasso regression
  ## For internals:
  #   data("Boston", package = "MASS")
  #   X <- model.matrix(medv~., Boston)
  #   y <- Boston$medv
  
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
                                search = "random", # NEEDS CHECKING: what is this search method?
                                selectionFunction = "best",
                                verboseIter = FALSE)
  
  # CV
  el_cv <- train(y ~., 
                 data = cbind(y, X), 
                 method = "glmnet",
                 family = glmfam, 
                 type.multinomial = "grouped",
                 trControl = train_control,
                 preProcess = c("center", "scale"),
                 tuneLength = 25
  )
  
  # Train model
  regu.mod <- glmnet(X, y,
                     family = glmfam,
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
  
  s2hat_mj <- mean((predict(model, X_tr) - y_tr)**2) 
    # according to paper this is the estimate
  X_te     <- model.matrix(rep(1, nrow(X_te)) ~., data.frame(X_te))[, -1]
    # create a prediction matrix (for possible dummy coded needed),
    # -1 no need for intercept
  yhat_te  <- rnorm(n = nrow(X_te),
                    mean = predict(model, X_te), 
                    sd = sqrt(s2hat_mj))
  z.m_j_mis <- yhat_te
  
  return(z.m_j_mis)
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
  py1       <- predict(model, X_te, type = "response") 
  z.m_j_mis <- rbinom(nrow(py1), 1, py1)
  # Clumsy way of fixing the levels
  z.m_j_mis <- factor(z.m_j_mis) # return to original labels
  levels(z.m_j_mis) <- levels(y_tr)
  
  return(z.m_j_mis)
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
  #   createDataPartition(p = 0.8, list = FALSE)
  # train  <- Boston[train_ind, ]
  # test <- Boston[-train_ind, ]
  # 
  # X_tr <- model.matrix(medv~., train)[,-1]
  # X_te <- model.matrix(medv~., test)[,-1]
  # y_tr <- train$medv
  # y_te <- test$medv
  #   cv <- cv.glmnet(X_tr, y_tr, alpha = 1) # cross validate lambda value
  # model <- glmnet(X_tr, y_tr, alpha = 1, lambda = cv$lambda.min)
  
  ## Body ##
  
  lasso.coef <- as.matrix(coef(model))
  S_hat <- row.names(lasso.coef)[lasso.coef != 0][-1] # estimated active set
  
  # Estimate Imputation model on observed z_j
  lm_fit <- lm(y_tr ~ X_tr[, S_hat])
  
  theta_MLE <- coef(lm_fit)
  Sigma_MLE <- vcov(lm_fit)
  
  # Sample parameters for prediction/imputation
  theta_m_j <- MASS::mvrnorm(1, 
                             theta_MLE, 
                             Sigma_MLE + diag(ncol(Sigma_MLE)) * parms$k_IURR )
  sigma_m_j <- sigma(lm_fit)
    # This is not correct. I'm using the OLS estiamtes with no variation in
    # sigma. Later fix this!
  
  # Get imputations
  y_imp <- rnorm(n = nrow(X_te),
                 mean = (model.matrix( ~ X_te[, S_hat]) %*% theta_m_j),
                 sd = sigma_m_j)
  
  return(y_imp)
}

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

get_50_best <- function(Xy_mis, S){
  ## Description
  # Given a dataset with missing univariate values, and the true active set
  # of the imputation model, it returns an index for the columns of the
  # dataset that selects the imputation target variable, the active set
  # and the first 50-length(active set) predictors for highest correlation
  # with the target variable.
  ## For internals
  var50_must <- c(S+1, ncol(Xy_mis))
  cortgt <- cor(Xy_mis[, -var50_must], use = "pairwise.complete.obs")[, "z1"]
  var50_opt <- sort(abs(cortgt), decreasing = TRUE)
  var50_opt <- var50_opt[2 : (50-length(S))]
  var50_opt <- which(colnames(Xy_mis) %in% names(var50_opt))
  var50 <- c(1, var50_must, var50_opt)
  
  return(var50)
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
