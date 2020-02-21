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
  a0 <- 1
  a <- rep(1, length(S)) * stnr
  z1 <- a0 * Z_m1[, S] %*% a
  
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

process_4regu <- function(Z, Zm, j_th, parms){
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

rr_est_lasso <- function(X, y, parms){
  ## Description:
  # Given any dv y (e.g. variable to be imputed), and its corresponding X 
  # values (observed, or imputed), fits a lasso regression
  ## For internals:
  #   data("Boston", package = "MASS")
  #   X <- model.matrix(medv~., Boston)
  #   y <- Boston$medv
  
  glmfam   <- detect_family(y) # regularized family type
  
  cv_lasso <- cv.glmnet(X, y,
                        family = glmfam,
                        type.multinomial = "grouped",
                        alpha = 1)
  
  regu.mod <- glmnet(X, y, 
                     family = glmfam,
                     type.multinomial = "grouped", # if glmfam is multinomial, otherwise does not apply
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

make_imp_gaus <- function(model, X_tr, y_tr, X_te, parms){
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

make_imp_dich <- function(model, X_tr, y_tr, X_te, parms){
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

make_imp_multi <- function(model, X_tr, y_tr, X_te, parms){
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

get_pool_est <- function(multi_dt, parms){
  ## Description:
  # Given a list of imputed datasets under the same imputation model
  # it returns the pooled estiamtes and CI of the regression coefs
  
  ## For internals
  # fake_data <- nhanes
  # colnames(fake_data) <- c("y", "z1", "z2", "z3")
  # imp <- mice(fake_data, print = FALSE, m = 10, seed = 24415) 
  # imp_dats <- mice::complete(imp, "all")
  # multi_dt = imp_dats
  
  ## Fit model
  models <- lapply(X = multi_dt,
                   FUN = function(x) lm(parms$formula, data = x))
  
  ## Extract info
  summa_models <- lapply(X = models,
                         FUN = function(x) summary(x))
  
  coefs <- t(sapply(X = summa_models,
                  FUN = function(x) coef(x)[, "Estimate"]))
  
  all_vcov <- lapply(X = summa_models,
                  FUN = function(x) vcov(x))
  
  ## Preliminary computations
  m <- length(multi_dt)
  
  Q_bar <- colMeans(coefs[1:(nrow(coefs)/2), ])

  U_bar <- diag(Reduce('+', all_vcov) / m)
  
  B <- diag( 1 / (m-1) * t(coefs - Q_bar) %*% (coefs - Q_bar) )
  
  T <- U_bar + B + B/m
  
  lambda <- (B+B/m)/T
  
  riv <- (B+B/m)/U_bar
  
  # degrees of freedom
  nu_old <- (m-1)*(1+1/riv**2)
  nu_com <- parms$n - parms$k
  nu_obs <- ( (nu_com+1) / (nu_com+3) ) * nu_com*(1-lambda)
  nu <- nu_old*nu_obs / (nu_old+nu_obs)
  
  ## CI computation
  t_nu <- qt((1-parms$alphaCI)/2, 
             df = nu, 
             lower.tail = FALSE)
  
  CI <- data.frame(lwr = Q_bar - t_nu * sqrt(T), 
                   upr = Q_bar + t_nu * sqrt(T))
  
  ## Output
  return(list(est = Q_bar,
              CI = CI))
}
