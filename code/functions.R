### Title:    helper functions
### Porject:  Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-05-19

# generic functions -------------------------------------------------------

update_report <- function(method_name = "Method ...", 
                          rep_count = 1, 
                          parms, 
                          cnd,
                          perform = TRUE){
  ## Exanple inputs
  # method_name = "Method ..."
  # rep_count = 1
  # cnd = data.frame(p = 500, pm = .3)
  # perform = TRUE
  # Updates the report .txt file with the count and name of accomplished task
  if(perform == TRUE){
    cat(paste0(Sys.time(), 
               " | Rep = ", rep_count, 
               ", ", paste0(names(cnd), " = ", cnd, collapse = ", "),
               " | ", method_name, 
               " done",
               "\n"),
        file = paste0(parms$outDir, parms$report_file_name),
        append = TRUE)
  }
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
  #     @Z0 <- mice::boys
  #     @missVarInfo <- missing_type(Z)
  # Output: a dataset with cotninuous variables imputed with mean value, and categorical with mode category
  # Used in: MICERandomForest
  # Notes: integer and numeric variables are considered (inaccurately) both as continuous
  ## Input examples from simulation
  # Z0 <- Z
  # missVarInfo <- missing_type(Z)
  
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

# Debugging function
# For use see: https://stackoverflow.com/questions/40629715/how-to-show-error-location-in-trycatch/40674718

withErrorTracing = function(expr, silentSuccess=FALSE) {
  hasFailed = FALSE
  messages = list()
  warnings = list()
  
  errorTracer = function(obj) {
    
    # Storing the call stack 
    calls = sys.calls()
    calls = calls[1:length(calls)-1]
    # Keeping the calls only
    trace = limitedLabels(c(calls, attr(obj, "calls")))
    
    # Printing the 2nd and 3rd traces that contain the line where the error occured
    # This is the part you might want to edit to suit your needs
    print(paste0("Error occuring: ", trace[length(trace):1][2:3]))
    
    # Muffle any redundant output of the same message
    optionalRestart = function(r) { res = findRestart(r); if (!is.null(res)) invokeRestart(res) }
    optionalRestart("muffleMessage")
    optionalRestart("muffleWarning")
  }
  
  vexpr = withCallingHandlers(withVisible(expr),  error=errorTracer)
  if (silentSuccess && !hasFailed) {
    cat(paste(warnings, collapse=""))
  }
  if (vexpr$visible) vexpr$value else invisible(vexpr$value)
}

# Estimation --------------------------------------------------------------

rr_est_ridge <- function(X, y, parms, fam="gaussian"){
  ## Description:
  # Given any dv y (e.g. variable to be imputed), and its corresponding X 
  # values (observed, or imputed), fits a lasso regression
  ## For internals:
  # data("Boston", package = "MASS")
  # X <- model.matrix(medv~., Boston)
  # y <- Boston$medv
  # fam <- "gaussian"
  ## Internals from simualtion
  # X = X_obs_bs
  # y = y_obs_bs
  # parms = parms
  # fam = glmfam
  
  cv_lasso <- cv.glmnet(X, y,
                        family = fam,
                        nfolds = 10,
                        alpha = 0)
  
  regu.mod <- glmnet(X, y, 
                     family = fam,
                     alpha = 1, 
                     lambda = cv_lasso$lambda.min)

  return(regu.mod)
}

rr_est_lasso <- function(X, y, parms, fam="gaussian"){
  ## Description:
  # Given any dv y (e.g. variable to be imputed), and its corresponding X 
  # values (observed, or imputed), fits a lasso regression
  ## For internals:
    # data("Boston", package = "MASS")
    # X <- model.matrix(medv~., Boston)
    # y <- Boston$medv
    # fam <- "gaussian"
  ## Internals from simualtion
  # X = X_obs_bs
  # y = y_obs_bs
  # parms = parms
  # fam = glmfam
  
  cv_lasso <- cv.glmnet(X, y,
                        family = fam,
                        nfolds = 10,
                        alpha = 1)
  
  regu.mod <- glmnet(X, y, 
                     family = fam,
                     alpha = 1, 
                     lambda = cv_lasso$lambda.min)
  
  return(regu.mod)
}

rr_est_elanet <- function(X, y, parms, fam = "gaussian"){
  # Source for cross validation strategy:
  # https://daviddalpiaz.github.io/r4sl/elastic-net.html
  ## Description:
  # Given any dv y (e.g. variable to be imputed), and its corresponding X 
  # values (observed, or imputed), fits a lasso regression
  ## For internals:
  # data("Boston", package = "MASS")
  # X <- model.matrix(medv~., Boston)
  # y <- Boston$medv
  # fam <- "gaussian"
  
  ## Internals from simualtion
  # X = X_obs_bs
  # y = y_obs_bs
  # parms = parms
  # fam = glmfam
  
  # Set training control
  train_control <- trainControl(method = "cv",
                                number = 10, # 10-fold cross validation 
                                selectionFunction = "best",
                                verboseIter = FALSE)
  
  # CV
  el_cv <- train(y ~ .,
                 data = data.frame(y, X), 
                 method = "glmnet",
                 family = fam, 
                 type.multinomial = "grouped",
                 trControl = train_control,
                 preProcess = c("center", "scale"),
                 tuneLength = 10 # values tried for both alpha and lambda
  )
  
  # Train model
  regu.mod <- glmnet(X, y,
                     family = fam,
                     type.multinomial = "grouped", # if glmfam is multinomial, otherwise does not apply
                     alpha = el_cv$bestTune$alpha,
                     lambda = el_cv$bestTune$lambda)
  
  return(regu.mod)
}

CFA_mod_wirte <- function(dat, lv_number, parms){
  # Define CFA model
  # lv_number = 2
  item_names <- colnames(dat)[1:(lv_number*parms$n_it)]
  lv_names <- paste0("lv", 1:lv_number)
  
  lv_models <- sapply(1:lv_number, function(i){
    paste0(lv_names[i], 
           " =~ ",
           paste0(item_names[((0:lv_number)[i]*parms$n_it+1):
                          ((0:lv_number)[i+1]*parms$n_it)], 
                  collapse = " + ")
    )
    
  })
  
  CFA_model <- paste(lv_models, collapse = "\n")
  
  return(CFA_model)
}

# SAT_mod_write <- function(parms){
#   var_names <- paste0("z", parms$z_m_id)
#   # Means
#   head_means <- "# Intercepts\n"
#   all_means <- paste0(var_names, " ~ ", "1")
#   all_means <- paste(all_means, collapse = "\n")
#   
#   # Variances
#   head_vars <- "# Variances \n"
#   all_vars <- paste0(var_names, " ~~ ", var_names)
#   all_vars <- paste(all_vars, collapse = "\n")
#   
#   # Coivariances
#   head_covs <- "# Covariances \n"
#   all_covs <- combn(var_names, 2)
#   all_covs <- apply(all_covs, 2, paste0, collapse = " ~~ ")
#   all_covs <- paste(all_covs, collapse = "\n")
#   
#   # Put together
#   SAT_mod <- paste(head_means,all_means,
#                    head_vars, all_vars,
#                    head_covs, all_covs
#   )
#   
#   return(SAT_mod)
# }

SAT_mod_write <- function(var_id){
  # var_id <- paste0("z", 1:10)
  # var_id <- colnames(SC_dt_sn$GS)[1:parms$sc_n]
  
  # Means
  head_means <- "# Means\n"
  all_means <- paste0(var_id, " ~ ", "1")
  all_means <- paste(all_means, collapse = "\n")
  
  # Variances
  head_vars <- "# Variances \n"
  all_vars <- paste0(var_id, " ~~ ", var_id)
  all_vars <- paste(all_vars, collapse = "\n")
  
  # Coivariances
  head_covs <- "# Covariances \n"
  all_covs <- combn(var_id, 2)
  all_covs <- apply(all_covs, 2, paste0, collapse = " ~~ ")
  all_covs <- paste(all_covs, collapse = "\n")
  
  # Put together
  SAT_mod <- paste(head_means,all_means,
                   head_vars, all_vars,
                   head_covs, all_covs
  )
  
  return(SAT_mod)
}


# Imputation --------------------------------------------------------------

imp_gaus_DURR <- function(model, X_tr, y_tr, X_te, parms){
  ## Description ##
  
  # Given a fitted imputation model it returns the imputations for
  # the missing values on a normally distributed imputation dv
  # (based on DengEtAl 2016, appendix)
  
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
  # the missing values on dichotomous. (based on DengEtAl 2016, appendix)
  
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
  
  ## Internals from simulation
  # model = regu.mod
  # X_tr = X_obs_bs
  # y_tr = y_obs_bs
  # X_te = X_mis
  
  ## Body ##
  
  X_te <- model.matrix( ~ X_te)[, -1]
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
  # according to IURR method (based on DengEtAl 2016, appendix)
  
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
  
  ## Inputs from simulation
  # model = regu.mod
  # X_tr = X_obs
  # y_tr = y_obs
  # X_te = X_mis
  # y_te = y_mis
  
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

imp_dich_IURR <- function(model, X_tr, y_tr, X_te, parms){
  ## Description ##
  
  # Given a regularized model it returns the imputations for
  # the missing values on bernulli distributed imputation dv
  # according to IURR method (based on DengEtAl 2016, appendix)
  
  ## For internals ##
  
  # data("Boston", package = "MASS")
  # train_ind <- Boston$medv %>%
  #   caret::createDataPartition(p = 0.8, list = FALSE)
  # train  <- Boston[train_ind, ]
  # test <- Boston[-train_ind, ]
  # 
  # X_tr <- model.matrix(chas~., train)[,-1]
  # X_te <- model.matrix(chas~., test)[,-1]
  # y_tr <- train$chas
  # y_te <- test$chas
  #   cv <- cv.glmnet(X_tr, y_tr, alpha = 1) # cross validate lambda value
  # model <- glmnet(X_tr, y_tr, alpha = 1, lambda = cv$lambda.min)
  # fam <- "binomial"
  
  ## Inputs from simulation
  # model = regu.mod
  # X_tr = X_obs
  # y_tr = y_obs
  # X_te = X_mis
  # y_te = y_mis
  
  ## Body ##
  
  # Select predictors based on rr
  rr_coef <- as.matrix(coef(model)) # regularized regression coefs
  rr_coef_no0 <- row.names(rr_coef)[rr_coef != 0]
  AS <- rr_coef_no0[-1] # predictors active set
  
  # 1. Obtain estiamtes w/ standard inference procedure (IWLS or ML)
  if(identical(rr_coef_no0, "(Intercept)")){
    glm_fit <- glm(y_tr ~ 1, family = "binomial")
  } else {
    glm_fit <- glm(y_tr ~ X_tr[, AS], family = "binomial")
  }

  # 2. obtain estimates
  theta <- coef(glm_fit)
  Sigma <- vcov(glm_fit)
  
  # 3. Sample parameters for posterior predictive distribution
  # (approximate distirbution of parameters)
  pdraws_par <- MASS::mvrnorm(1,
                              mu = theta, 
                              Sigma = Sigma)
  
  # 4. Sample approxiamtion posterior predictive distribution
  if(identical(rr_coef_no0, "(Intercept)")){
    logit <- b_ppd <- pdraws_par # only intercept
    prob <- exp(logit) / (1+exp(logit))
    y_imp <- rbinom(nrow(X_te), 1, prob)
  } else {
    X_ppd <- model.matrix( ~ X_te[, AS]) # X for posterior pred dist
    b_ppd <- pdraws_par # betas for posterior pred dist
    logit <- X_ppd %*% b_ppd
    prob <- exp(logit) / (1+exp(logit))
    y_imp <- rbinom(nrow(X_te), 1, prob)
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

.fmi_compute <- function(fits){
  ## Description
  # Given a list of fits on multiply imputed datasets
  # it returns the fmi for each estiamted paramter
  
  ## For internals
  # fits = semR_fit_mi[["bridge"]]
  
  ## Body
  ## Coef estimates
  m <- length(fits)

  coefs <- sapply(X = fits,
                  FUN = function(x) coef(x))
  Q_bar <- rowMeans(coefs)
  
  ## Variances (squared standard error for each parameter)
  all_vcov <- lapply(X = fits,
                     FUN = function(x) vcov(x))
  U_bar <- diag(Reduce('+', all_vcov) / m)
  
  B <- diag(1 / (m-1) * (coefs - Q_bar) %*% t(coefs - Q_bar))
  
  T_var <- U_bar + B + B/m
  
  # FMI vector
  FMI <- round(.fmi(m = m, b = B, t = T_var), 3)
  
  return(FMI)
}

# Estimate regression coefficeints

# fit_lm_models <- function(multi_dt, vrbs){
#   ## Description:
#   # Given a list of complete datasets it fits a linear model
#   # to obtain standardized regression coefficients (all vairables 
#   # are centered and standardized)
#   ## Example internals
#   # multi_dt <- imp_DURR_la$dats
#   # vrbs <- parms$lm_model
#   if(!is.null(multi_dt)){
#   mod <- paste0(vrbs[1], 
#             " ~ - 1 + ", 
#             paste0(vrbs[-1], collapse = " + ")
#   )
#     models <- lapply(X = multi_dt,
#                      FUN = function(x) lm(mod, 
#                                           data = as.data.frame( scale(x) )
#                                           )
#                      )
#   } else {models = NULL}
#   return(models)
# }

fit_lm_models <- function(multi_dt, mod){
  ## Description:
  # Given a list of complete datasets it fits a linear model
  # to obtain standardized regression coefficients (all vairables 
  # are centered and standardized)
  ## Example internals
  # multi_dt <- imp_DURR_la$dats
  # vrbs <- parms$lm_model
  if(!is.null(multi_dt)){
    models <- lapply(multi_dt,
                     function(x) lm(mod, 
                                    data = as.data.frame( x )
                     )
    )
  } else {models = NULL}
  return(models)
}

# MLE Estiamtes of means, variances, covariances

fit_sat_model <- function(multi_dt){
  # Given a list of complete datasets it fits a model described
  # in mod
  ## Example input ##
  # multi_dt <- list(imp_PCA$dats,
  #                  imp_MICE_TR$dats)[[1]]
  ## Body ##
  if(!is.null(multi_dt)){
    models <- lapply(X = multi_dt,
                     FUN = function(x) {
                       tryCatch({
                         # Obtain MLE estiamtes
                         sem(parms$lav_model, 
                             data = x, 
                             likelihood = "wishart")},
                         # If there is a fitting error, report it
                         error = function(report) {
                           err <- paste0("Original Error: ", report)
                           return(err)
                         },
                         # if there is a warning error, report it
                         warning = function(report) {
                           err <- paste0("Original Warning: ", report)
                           return(err)
                         })
                     })
    
    # Keep only models that converged
    fits_indx <- as.vector(which(!sapply(models, is.character)))
    models <- models[fits_indx]
    
  } else {models = NULL}
  return(models)
}

sem_EST <- function(fits){
  ## Description
  # Given a list of fits for different single imputation apoprahces
  # it returns the estiamtes of the parameters for each
  ## For internals
  # fits = list(fit_mf, fit_gs, fit_cc)
  # summa_models <- lapply(X = fits,
  #                        FUN = function(x) parameterEstimates(x))
  # 
  # coefs <- sapply(X = summa_models,
  #                 FUN = function(x) x[, c("est")])
  coefs <- sapply(X = fits,
                  FUN = function(x) coef(x))
  return(coefs)
}

sem_CI <- function(fits){
  ## Description
  # Given a list of fits for different single imputation apoprahces
  # it returns the CI of the estiamtes of the parameters for each
  ## For internals
  # fits = list(fit_mf, fit_gs, fit_cc)
  
  row_indx <- !is.na(parameterEstimates(fits[[1]])[, "z"])
    # In the CFA model, either the variances of lv or factor loadings
    # are fixed to some value. Hence, they would not be paramters.
    # However, while the vcov function keeps this into account, the
    # computation of Q_bar that starts from the output of paramterEstimates
    # does not. So I need an index that identifies what rows of
    # paramterEstimates output are actual paramters.
  CI_list <- lapply(X = fits,
                    FUN = function(x) parameterEstimates(x)[row_indx, 
                                                            c("ci.lower", 
                                                              "ci.upper")]
  )
  
  CI_mtx <- sapply(CI_list, function(x){
    c(x[,1], x[,2])
  })

  return(CI_mtx)
}

lm_EST <- function(fits){
  ## Description
  # Given a list of fits for different single imputation apoprahces
  # it returns the estiamtes of the parameters for each
  ## For internals
  # fits = lm_sndt
  
  coefs <- sapply(X = fits,
                         FUN = function(x) coef(x))
  
  return(coefs)
}

lm_CI <- function(fits){
  ## Description
  # Given a list of fits for different single imputation apoprahces
  # it returns the CI of the estiamtes of the parameters for each
  ## For internals
  # fits = lm_sndt
  CI_list <- lapply(X = fits, confint)
                    
  CI_mtx <- sapply(CI_list, function(x){
    c(x[,1], x[,2])
  })
  
  return(CI_mtx)
}
  
sem_pool_EST_f <- function(fits){
  ## Description
  # Given a list of fits from different imputed datasets under the
  # same imputation model it returns the pooled estiamtes of the
  # parameters of interest.
  ## For internals
  # fits = SAT_fit_sc_mi[[1]] # fits_md[[4]]
  ## Pool coefs
  coefs <- sapply(X = fits,
                  FUN = function(x) coef(x))
  Q_bar <- rowMeans(coefs)
    # 1 column per imputed dataset
  return(Q_bar)
}

sem_pool_CI_f <- function(fits){
  ## Description
  # Given a list of imputed datasets under the same imputation model
  # it returns the pooled CIs of the regression coefs
  ## Note on storing convention: instead of storing low and up bound
  # in different columns I stored in one single column. The first half
  # is lwr bound, the bottom part is upper bound. This makes it easier
  # to store. 
  
  ## For internals
  # fits = SAT_fit_raw_mi[[1]]
  # fits = sem_fits[lapply(sem_fits, length) != 0][[1]]
  
  ## Body
  ## Coef estimates
  m <- length(fits)
  # use only dataset for which model can be fit
  row_indx <- !is.na(parameterEstimates(fits[[1]])[, "z"])
    # In the CFA model, either the variances of lv or factor loadings
    # are fixed to some value. Hence, they would not be paramters.
    # However, while the vcov function keeps this into account, the
    # computation of Q_bar that starts from the output of paramterEstimates
    # does not. So I need an index that identifies what rows of 
    # paramterEstimates output are actual paramters.
  
  # summa_models <- lapply(X = fits,
  #                        FUN = function(x) parameterEstimates(x)[row_indx,])
  # coefs <- sapply(X = summa_models,
  #                 FUN = function(x) x[, c("est")])
  coefs <- sapply(X = fits,
                  FUN = function(x) coef(x))
  Q_bar <- rowMeans(coefs)
  
  ## Variances (squared standard error for each parameter)
  all_vcov <- lapply(X = fits,
                     FUN = function(x) vcov(x))
  U_bar <- diag(Reduce('+', all_vcov) / m)
  
  B <- diag(1 / (m-1) * (coefs - Q_bar) %*% t(coefs - Q_bar))

  T_var <- U_bar + B + B/m
  
  ## Degrees of freedom
  nu_com <- parms$n - nrow(coefs) # n - k where k number of paramteres estimated
  nu <- .miDf(length(fits), b = B, t = T_var, nu_com)
 
  ## CI computation
  t_nu <- qt(1 - (1-parms$alphaCI)/2, 
             df = nu)
  
  CI <- c(lwr = Q_bar - t_nu * sqrt(T_var), 
                   upr = Q_bar + t_nu * sqrt(T_var))
  
  return(CI = CI)
}

## LM pooling

lm_pool_EST_f <- function(fits){
  ## Description
  # Given a list of imputed datasets under the same imputation model
  # it returns the pooled estiamtes of the regression coefs
  ## For internals
  # fits <- lm_fits[[2]]
  
  # Extract estiamtes from the fitted models
  summa_models <- lapply(X = fits,
                         FUN = function(x) summary(x))
  coefs <- t(sapply(X = summa_models,
                    FUN = function(x) coef(x)[, "Estimate"]))
  
  Q_bar <- colMeans(coefs)
  
  return(Q_bar)
  
}

lm_pool_CI_f <- function(fits){
  ## Description
  # Given a list of imputed datasets under the same imputation model
  # it returns the pooled CIs of the regression coefs
  
  ## Example internals
  # fits <- lm_fits[[1]]
  
  # Do we have fits for a given imputation method?
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
  
  CI <- c(lwr = Q_bar - t_nu * sqrt(T), 
          upr = Q_bar + t_nu * sqrt(T))
  
  return(CI = CI)
}

onetree <- function(xobs, xmis, yobs, s) {
  ## Ripped off mice package mice.impute.rf. Fits one free for the 
  ## random forest. Used in your own random forest version
  fit <- randomForest::randomForest(x = xobs, y = yobs,
                                    ntree = 1)
  leafnr <- predict(object = fit, newdata = xobs, nodes = TRUE)
  nodes <- predict(object = fit, newdata = xmis, nodes = TRUE)
  donor <- lapply(nodes, function(s) yobs[leafnr == s])
  return(donor)
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

## NEW VERSIONS 
fit_sem <- function(multi_dt, model, std.lv=FALSE){
  # Given a list of complete datasets it fits a model described
  # in mod
  ## Example input ##
  # multi_dt <- imp_DURR_la$dats
  # multi_dt <- SC_dt_sn
  # model <- SAT_mod_write(colnames(SC_dt_sn$GS)[1:parms$sc_n],
  #                        parms, score = TRUE)
  ## Body ##
  if(!is.null(multi_dt)){
    fits <- lapply(X = multi_dt,
                     FUN = function(x) {
                       tryCatch({
                         # Fit SEM model
                         sem(model = model, 
                             data = x, 
                             likelihood = "wishart",
                             std.lv = std.lv)},
                         # If there is a fitting error, report it
                         error = function(report) {
                           err <- paste0("Original Error: ", report)
                           return(err)
                         },
                         # if there is a warning error, report it
                         warning = function(report) {
                           err <- paste0("Original Warning: ", report)
                           return(err)
                         })
                     })
    # Keep only models that converged (no fitting errors)
    fits_indx <- as.vector(which(!sapply(fits, is.character)))
    fits <- fits[fits_indx]
    
  } else {fits = NULL}
  return(fits)
}

fit_lm <- function(multi_dt, model){
  ## Description:
  # Given a list of complete datasets it fits a linear model
  # to obtain standardized regression coefficients (all vairables 
  # are centered and standardized)
  ## Example internals
  # multi_dt = SC_dt_mi$DURR_la
  # model = LM_formula
  if(!is.null(multi_dt)){
    fits <- lapply(X = multi_dt,
                     FUN = function(x) lm(model, 
                                          data = x)
    )
  } else {fits = NULL}
  return(fits)
}

# Create Scores for MI datasets
scorify <- function(dat_in, cond, parms){
  ## Makes raw data into univariate scores corresponding to theoretical 
  ## constructs.
  # Example inputs
  # dat_in <- Xy # single dataset
  
  dat_out <- data.frame(matrix(nrow = nrow(dat_in), ncol = cond$lv))
    colnames(dat_out) <- paste0("sc", 1:cond$lv)
  
  for (i in 1:cond$lv) {
    item_idx <- c((0:cond$lv)[i]*parms$n_it+1):((0:cond$lv)[i+1]*parms$n_it)
    dat_out[, i] <- rowMeans(dat_in[, item_idx])
  }
  
  return(dat_out)
    
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

mean_traceplot <- function(out, 
                           dat = 1, # which data repetition should I show?
                           method = "blasso", # same name as in parms
                           y_range = c(-10, 20),
                           iters = 1:5){
  ## Internals
  # out <- out_cnv
  # iters = 200:300

  ## Description
  # It prints the traceplots for the mean imputed values in each iteration
  # in different chains, by variable, one dataset, one imputation method
  
  # Display in same pane
  par(mfrow = c(3, ceiling(out$parms$zm_n/3)))
  # Plot
  # Are imputations of mids class?
  if(class(out[[dat]][[1]]$imp_values[[method]]) == "mids"){
    plot(out[[dat]][[1]]$imp_values[[method]])
  } else {
    for (v in 1:length(out$parms$z_m_id)) {
      # CHAIN 1
      # Mean imputed value across individuals in each iteration
      mean_imp <- rowMeans(out[[dat]][[1]]$imp_values[[method]][[1]][[v]][iters, ])
      
      plot(iters, mean_imp, type = "l",
         main = method,
         ylim = y_range,
         ylab = paste0("z", v), xlab = "Iteration")
      # plot(seq(out$parms$iters)[iters], mean_imp, type = "l",
      #      main = method,
      #      ylim = y_range,
      #      ylab = paste0("z", v), xlab = "Iteration")
      
      # CHAIN 2 to m 
      for (i in 2:(out$parms$chains)) {
        mean_imp <- rowMeans(out[[dat]][[1]]$imp_values[[method]][[i]][[v]][iters, ])
        lines(iters, mean_imp)
      }
    }
  }
}
  
res_sem_time <- function(out, condition = 1){
  # Sem Model
  select_cond <- names(out[[1]])[condition]
  
  # Time
  res_time <- NULL
  for (i in 1:out$parms$dt_rep) {
    res_time <- rbind(res_time, out[[i]][[select_cond]]$run_time_min)
  }
  return( round(colMeans(res_time), 3) )
}

res_sum <- function(out, model, condition = 1, bias_sd = FALSE){
  # model = "semS" # the first part of the name of a result object stored
                   # in out 
  # model = "lm"
  # condition = 1
  
  ## Prep ##
  est <- paste0(model, "_EST")
  ci <- paste0(model, "_CI")
  select_cond <- names(out[[1]])[condition]
  
  ## Step 2. Bias ##
  avg <- sapply(out$parms$methods, function(m){
    store <- NULL
    for (i in 1:out$parms$dt_rep) {
      succ_method <- colnames(out[[i]][[select_cond]][[est]])
      store <- cbind(store, 
                     out[[i]][[select_cond]][[est]][, 
                                                     succ_method %in% m])
    }
    c(rowMeans(store, na.rm = TRUE), rep = ncol(store)) # MCMC statistics 
  })
  # Store Objects
  validReps  <- avg["rep", ] # number of successes
  avg        <- avg[-which(rownames(avg) == "rep"), ]
  psd_tr_vec <- avg[, "GS"] # pseudo true values
  
  if(out$parms$exp == 1 & nchar(rownames(avg)[1])<1){
    # Fixes a problem of old results from privous runs of eperiment 1
    # in the future I will delte these part as exp 1 results wiil be
    # uniform with rest
    fit <- lavaan::sem(out$parms$lav_model,
                       data = out[[1]][[select_cond]]$dat_full,
                       likelihood = "wishart")
    rownames(avg) <- apply(parameterEstimates(fit)[,1:3], 
                           1, 
                           paste0, 
                           collapse = "")
  }
  
  # Raw bias
  bias <- avg - psd_tr_vec
  
  # Bias as percentage of true value
  bias_per <- cbind(ref = round(psd_tr_vec, 3),
                    round(
                      abs(bias)/psd_tr_vec*100, 
                      0)
  )
  
  meths <- out$parms$methods[-which(out$parms$methods == "GS")]
  # Bias Mean Standardized
  if(bias_sd == TRUE){
    means_indx <- grep("~1", rownames(avg))
    vars_indx <- means_indx + length(means_indx)
    bias_sd <- round(
      (avg[means_indx, ] - psd_tr_vec[means_indx])/
        sqrt(psd_tr_vec[vars_indx]),
      3)
  } else {
    bias_sd <- NULL
  }
  
  ## Step 3. CI Coverange ##
  # storing threshold
  str_thrs <- nrow(out[[1]][[select_cond]][[ci]])/2
  
  # Confidence Interval Coverage
  CIC <- sapply(out$parms$methods, function(m){
    store <- NULL
    for (i in 1:out$parms$dt_rep) {
      succ_method <- colnames(out[[i]][[select_cond]][[est]])
      col_indx <- succ_method %in% m
      cond_est <- out[[i]][[select_cond]][[est]]
      cond_CI  <- out[[i]][[select_cond]][[ci]]
      ci_low   <- cond_CI[1:str_thrs, ]
      ci_hig   <- cond_CI[-(1:str_thrs), ]
      
      store <- cbind(store, 
                     ci_low[, col_indx] < psd_tr_vec &
                       psd_tr_vec < ci_hig[, col_indx]
      )
    }
    rowMeans(store, na.rm = TRUE) # MCMC statistics 
  })
  rownames(CIC) <- rownames(bias)
  
  # Output
  results <- list(cond      = select_cond,
                  MCMC_est  = round(cbind(ref=psd_tr_vec, avg), 3),
                  bias_raw  = round(cbind(ref=psd_tr_vec, bias), 3),
                  bias_per  = bias_per,
                  bias_sd   = bias_sd,
                  ci_cov    = round(CIC*100, 1),
                  validReps = validReps)
  return(results)
}

res_sem_sum <- function(out, condition = 1){
  ## Prep ##
  select_cond <- names(out[[1]])[condition]
  
  ## Step 1. Obtain Pseudo True Values ##
  full_dat_est <- matrix(NA, 
                         # Data repetitions
                         nrow = out$parms$dt_rep, 
                         # Parameters estiamtes
                         ncol = nrow(out[[1]][[select_cond]]$sem_EST))
  for (i in 1:out$parms$dt_rep) {
    full_dat_est[i, ] <- out[[i]][[select_cond]]$sem_EST[, which(out$parms$methods == "GS")]
  }
  
  psd_tr_vec <- colMeans(full_dat_est) # pseudo true values
  
  ## Step 2. Bias ##
  avg <- sapply(out$parms$methods, function(m){
    store <- NULL
    for (i in 1:out$parms$dt_rep) {
      succ_method <- colnames(out[[i]][[select_cond]]$sem_EST)
      store <- cbind(store, 
                     out[[i]][[select_cond]]$sem_EST[, 
                                                     succ_method %in% m])
    }
    c(rowMeans(store, na.rm = TRUE), rep = ncol(store)) # MCMC statistics 
  })
  
  # Store Objects
  validReps <- avg["rep", ] # number of successes
  avg <- avg[-which(rownames(avg) == "rep"), ]
  
  # Define names of parameters
  fit <- lavaan::sem(out$parms$lav_model,
                     data = out[[1]][[select_cond]]$dat_full,
                     likelihood = "wishart")
  rownames(avg) <- apply(parameterEstimates(fit)[,1:3], 
                         1, 
                         paste0, 
                         collapse = "")

  # Raw bias
  bias <- avg - psd_tr_vec

  # Bias as percentage of true value
  bias_per <- cbind(ref = round(psd_tr_vec, 3),
                    round(
                      abs(bias)/psd_tr_vec*100, 
                      0)
                    )
  
  ## Step 3. CI Coverange ##
  # storing threshold
  str_thrs <- nrow(out[[1]][[select_cond]]$sem_CI)/2
  
  # Confidence Interval Coverage
  CIC <- sapply(out$parms$methods, function(m){
    store <- NULL
    for (i in 1:out$parms$dt_rep) {
      succ_method <- colnames(out[[i]][[select_cond]]$sem_EST)
      col_indx <- succ_method %in% m
      cond_est <- out[[i]][[select_cond]]$sem_EST
      cond_CI  <- out[[i]][[select_cond]]$sem_CI
      ci_low   <- cond_CI[1:str_thrs, ]
      ci_hig   <- cond_CI[-(1:str_thrs), ]
      
      store <- cbind(store, 
                     ci_low[, col_indx] < psd_tr_vec &
                       psd_tr_vec < ci_hig[, col_indx]
                     )
    }
    rowMeans(store, na.rm = TRUE) # MCMC statistics 
  })
  rownames(CIC) <- rownames(bias)
  
  # Output
  results <- list(cond = select_cond,
                  MCMC_est = round(cbind(ref=psd_tr_vec, avg), 3),
                  bias_raw = round(cbind(ref=psd_tr_vec, bias), 3),
                  bias_per = bias_per,
                  ci_cov   = round(CIC*100, 1),
                  validReps = validReps)
  return(results)
}

res_lm_sum <- function(out, condition = 1){
  # Sem Model
  select_cond <- names(out[[1]])[condition]
  
  ## 1. Obtain Pseudo True Values ##
  
  full_dat_est <- matrix(NA, 
                         nrow = out$parms$dt_rep, 
                         ncol = nrow(out[[1]][[select_cond]]$lm_EST))
  for (i in 1:out$parms$dt_rep) {
    full_dat_est[i, ] <- out[[i]][[select_cond]]$lm_EST[, which(out$parms$methods == "GS")]
  }
  
  psd_tr_vec <- colMeans(full_dat_est) # pseudo true values
  
  ## Step 2. Compute averages of statistics (MCMC estiamtes) ##
  
  ## Step 2. Bias ##
  avg <- sapply(out$parms$methods, function(m){
    store <- NULL
    for (i in 1:out$parms$dt_rep) {
      succ_method <- colnames(out[[i]][[select_cond]]$lm_EST)
      store <- cbind(store, 
                     out[[i]][[select_cond]]$lm_EST[, 
                                                     succ_method %in% m])
    }
    c(rowMeans(store, na.rm = TRUE), rep = ncol(store)) # MCMC statistics 
  })
  
  # Store Objects
  validReps <- avg["rep", ] # number of successes
  avg <- avg[-which(rownames(avg) == "rep"), ]
  MCMC_est  <- round(cbind( ref = psd_tr_vec, avg), 3)
  
  # Raw bias
  bias_raw <- avg - psd_tr_vec
  
  # Bias as percentage of true value
  bias_per <- cbind(ref = round(psd_tr_vec, 3),
                    round(
                      abs(bias_raw)/psd_tr_vec*100, 
                      0)
  )
  
  ## Step 3: Obtain CI Coverages ##
  # storing threshold
  str_thrs <- nrow(out[[1]][[select_cond]]$lm_CI)/2
  
  # Confidence Interval Coverage
  CIC <- sapply(out$parms$methods, function(m){
    store <- NULL
    for (i in 1:out$parms$dt_rep) {
      succ_method <- colnames(out[[i]][[select_cond]]$lm_EST)
      col_indx <- succ_method %in% m
      cond_est <- out[[i]][[select_cond]]$lm_EST
      cond_CI  <- out[[i]][[select_cond]]$lm_CI
      ci_low   <- cond_CI[1:str_thrs, ]
      ci_hig   <- cond_CI[-(1:str_thrs), ]
      
      store <- cbind(store, 
                     ci_low[, col_indx] < psd_tr_vec &
                       psd_tr_vec < ci_hig[, col_indx]
      )
    }
    rowMeans(store, na.rm = TRUE) # MCMC statistics 
  })
  
  # Output
  results <- list(cond = select_cond,
                  MCMC_est = MCMC_est,
                  bias_raw = bias_raw,
                  bias_per = bias_per,
                  ci_cov = CIC*100,
                  validReps = validReps)
  return(results)
}

res_ed_est <- function(results, index = 1:2){
  # Internals
  # results = exp2_res$semR$cond1
  # index   = 1:2
  
  # Prepare objects for distance computation
  method_id <- colnames(results$MCMC_est) != "ref"
  ref <- results$MCMC_est[index, "ref"]
  MCMC_est <- results$MCMC_est[index, method_id]
  
  # Compute Euclidean distance
  out_dist <- sapply(as.data.frame(MCMC_est), 
                     function(x) dist(rbind(x, ref),
                                      method = "euclidean")
  )
  
  # Prepare and return output
  return( round(out_dist, 3) )
}

res_ed_ci <- function(results){
  # Internals
  # results <- semR_res[[1]]
  
  # Prepare objects for distance computation
  ref <- rep(95, nrow(results$ci_cov))
  
  # Compute Euclidean distance
  out_dist <- sapply(as.data.frame(results$ci_cov), 
                     function(x) dist(rbind(x, ref),
                                      method = "euclidean")
  )
  
  # Prepare and return output
  return( round(out_dist, 1) )
}

bridge_cv <- function(out, mods = NULL){
  # Returns a df with ridge penality selected for each condition
  # Compute Average FMI across all parameter estiamtes per ridge value
  
  # Arguments check
  if(is.null(mods)) mods = names(out[[1]][[1]]$fmi)
  
  # Body
  store_0 <- list()
  for (i in 1:nrow(out$conds)) {
    store_1 <- NULL
    for (dt in 1:out$parms$dt_rep) {
      store_1 <- cbind(store_1, unlist(out[[dt]][[i]]$fmi[mods]))
    }
    # store_0 <- cbind(store_0, rowMeans(store_1))
    store_0[[i]] <- rowMeans(store_1)
  }
  ridge_range <- length(unique(out$conds$ridge))
  names(store_0) <- rep(unique(out$conds$ridge), nrow(out$conds)/ridge_range)
  
  avg_fmi <- round(sapply(store_0, mean), 3)
  
  # Select ridge value with lowest average FMI
  i <- 1; j <- i+ridge_range-1
  ridge_s <- NULL
  for (r in 1:(length(avg_fmi)/ridge_range)) {
    ridge_s[r] <- as.numeric(names(which.min(avg_fmi[i:j])))
    i <- j+1
    j <- i+ridge_range-1
  }
  
  # Attach ridge value to specific condition
  col_indx <- colnames(out$conds) != "ridge" # exclude ridge column
  output <- data.frame(out$conds[!duplicated(out$conds[, col_indx]), 
                                 col_indx],
                       ridge = ridge_s)
  return(output)
}
