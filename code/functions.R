### Title:    helper functions
### Porject:  Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-05-19

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
  #     Z <- Z_mm
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
  # Z0 <- Z_mm
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
      Z0 <- Z0 %>% mutate_at(vars(c(missVarInfo$contVars, 
                                    missVarInfo$ordeVars)[j]),
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
  # X = X_obs
  # y = y_obs
  # parms = parms
  # fam = "gaussian"
  
  cv_lasso <- cv.glmnet(X, y,
                        family = fam,
                        nfolds = 10,
                        alpha = 1)
  
  return(cv_lasso)
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

exp4_fit_mod1 <- function(multi_dt){
  ## Description:
  # Given a list of complete datasets it fits a linear model
  # to obtain standardized regression coefficients (all vairables 
  # are centered and standardized)
  ## Example internals
  # multi_dt <- imp_DURR_la$dats
  # multi_dt <- si_data
  if(!is.null(multi_dt)){
    models <- lapply(multi_dt,
                     function(m){
                       # m <- multi_dt[["GS"]]
                       # Euthanasia
                       euth <- m$v156

                       # Country
                       country <- m$country
                       
                       # General Trust
                       trust.g <- m$v31
                       
                       # Confidence in Health care sys
                       trust.hs <- match(m$v126, 4:1)
                       
                       # Confidence in press
                       trust.pr <- m$v118
                       
                       # Confidence in state scale (need to create scale)
                       trust.s.var <- c("v121", "v120", "v127", "v131")
                       invCod <- sapply(trust.s.var, 
                                        function(x) match(m[, x], 4:1))
                       m[, trust.s.var] <- invCod
                       trust.s <- rowMeans(m[, trust.s.var]) # inverted
                       
                       # Age
                       age <- m$age
                       
                       # Education
                       edu <- m$v243_ISCED_1
                       
                       # Gender (1 = male)
                       sex <- m$v225
                       
                       # Religiousness
                       rel <- match(m$v6, 4:1)-1
                       
                       # Religious Denomination
                       denom <- m$v51v52_comb

                       # Fit mode
                       mod <- lm(euth ~
                                   # trust.g + 
                                   trust.hs + trust.pr + trust.s +
                                   edu + sex + age + rel + denom + country)
                       return(mod)
                     }
                     )
  } else {models = NULL}
  return(models)
}

exp4_fit_mod2 <- function(multi_dt){
  ## Description:
  # Given a list of complete datasets it fits a linear model
  # to obtain standardized regression coefficients (all vairables 
  # are centered and standardized) (ImmerzeelEtAl2016)
  ## Example internals
  # multi_dt <- imp_DURR_la$dats
  # multi_dt <- si_data
  if(!is.null(multi_dt)){
    models <- lapply(multi_dt,
                     function(m){
                       # Left / Right voting
                       lr <- m$v174_LR
                       
                       # Country
                       country <- m$country  

                       # Female
                       female <- relevel(m$v225, ref = "male")
                       
                       # Employment Status  
                       SES <- m$v246_egp
                       
                       # Native attitudes (mean of itmes)  
                       nativ <- c("v185", # jobs
                                  "v186", # crime
                                  "v187"  # strain on welfare
                       )
                       NatAt <- rowMeans(m[, nativ])
                       
                       # Authoritarian Attitudes
                       # Low and order attitudes
                       strongL <- match(m$v145, 4:1)
                       order <- forcats::fct_collapse(m$v110,
                                                      no = c("fighting rising prices", 
                                                             "more say in important government decisions", 
                                                             "protect freedom of speech"),
                                                      yes = c("maintaining order in nation")
                       )
                       
                       # Political Interest
                       polInterest <- match(m$v97, 4:1)
                       
                       # Political Action
                       action <- paste0("v", 98:101)
                       polAction <- rowMeans(m[, action])
                       polAction_r <- case_when(
                         polAction %in% 1 ~ 3, # have done all
                         polAction %in% seq(1.25, 2.75, by = .25) ~ 2,
                         polAction %in% 3 ~ 1 # would never do any
                       )
                       
                       # Covariates
                       age <- m$age
                       edu <- m$v243_ISCED_1
                       mat <- m$v234
                       
                       # urb <- m$v276_r
                       urb <- forcats::fct_collapse(m$v276_r,
                                                    less5e3 = c("under 5000"),
                                                    mid_5e3_2e4 = c("5000-20000"),
                                                    mid_2e4_1e5 = c("20000-100000"),
                                                    more1e5 = c("100000-500000", 
                                                                "500000 and more")
                       )
                       # Religiousness
                       rel <- match(m$v6, 4:1)-1
                       
                       # Denomination
                       denom <- m$v51v52_comb
                      
                       # Fit model
                       mod <- lm(lr ~ female + NatAt + strongL + polInterest +
                                   polAction_r + age + edu + rel + denom + country)
                       return(mod)
                     }
    )
  } else {models = NULL}
  return(models)
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
  # model <- cv.glmnet(X_tr, y_tr, alpha = 1) # cross validate lambda value
  # coef(model, s = "lambda.min")
  # model <- glmnet(X_tr, y_tr, alpha = 1, lambda = cv$lambda.min)
  # coef(model)
  
  ## Inputs from simulation
  # model = regu.mod
  # X_tr = X_obs
  # y_tr = y_obs
  # X_te = X_mis
  # y_te = y_mis
  
  ## Body ##
  
  # Select predictors based on rr
  rr_coef <- as.matrix(coef(model, s = "lambda.min")) # regularized regression coefs
  rr_coef_no0 <- row.names(rr_coef)[rr_coef != 0]
  AS <- rr_coef_no0[-1] # predictors active set
  
  # MLE estimate of model parameters
  # # ORIGINAL  
  # # 1. define starting values
  #   if(identical(rr_coef_no0, "(Intercept)")){
  #     lm_fit <- lm(y_tr ~ 1)
  #     X_mle <- model.matrix(y_tr ~ 1)
  #   } else {
  #     X_mle <- model.matrix(y_tr ~ X_tr[, AS])
  #     lm_fit <- lm(y_tr ~ X_tr[, AS])
  #   }
  #   startV <- c(coef(lm_fit), sigma(lm_fit))
    
  # UDPATE
    # 1. define starting values
    if(identical(rr_coef_no0, "(Intercept)")){
      lm_fit <- lm(y_tr ~ 1)
      X_mle <- model.matrix(y_tr ~ 1)
    } else {
      lm_fit <- lm(y_tr ~ X_tr[, AS])
      X_mle  <- model.matrix(y_tr ~ X_tr[, AS])
      colnames(X_mle) <- str_replace(colnames(X_mle), ".*]+", "")
      
      b.estimated <- coef(lm_fit)[!is.na(coef(lm_fit))]
      b.names <- str_replace(names(b.estimated), ".*]+", "")
      
      X_mle  <- X_mle[, colnames(X_mle) %in% b.names]
      # Fix NAs when coefficinet cannot be estiamted because variable
      # is near constant
    }
    
    # Option 1: Get rid of estimate
    startV <- c(coef(lm_fit)[!is.na(coef(lm_fit))], 
                # keep only coefficients for variables that were not kicked
                # out of the equation
                sigma(lm_fit))
    
    # # Option 2: Include but give custom 0 as initial value
    # startV <- c(coef(lm_fit), sigma(lm_fit))
    # startV[is.na(startV)] <- 0
    
  # 2. optimize loss function
    MLE_fit <- optim(startV, 
                     .lm_loss,
                     method = "BFGS",
                     hessian = T,
                     y = y_tr, X = X_mle)
    
  # 3. obtain estimates
    theta <- MLE_fit$par
    OI <- solve(MLE_fit$hessian) # parameters cov matrix

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
      colnames(X_ppd) <- str_replace(colnames(X_ppd), ".*]+", "")
      X_ppd  <- X_ppd[, colnames(X_ppd) %in% b.names]
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
  # source: https://rpubs.com/YaRrr/MLTutorial
  ## Example Inputs
  # X = X_mle
  # y = y_obs
  # theta = startV
  ## Body
  # prep
  k <- ncol(X)
  beta <- theta[1:k]
  sigma <- theta[k+1]
  
  # Check validity of sigma
  if(sigma < 0) { 
    # the optimization procedure to stay away from invalid parameter values.
    dev <- 1e7
  } else {
    # calculate (log) likelihood of each data point
    ll <- dnorm(y, mean = X %*% beta, sd = sigma, log = TRUE)
    # summarize into deviance score
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
                         # Obtain MLE estimates
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
  # fits <- m1_mi$blasso
  
  # Extract estiamtes from the fitted models
  coefs <- t(sapply(X = fits,
                    FUN = function(x) coef(x)))
  
  # Take mean
  Q_bar <- colMeans(coefs)
  
  return(Q_bar)
}

lm_pool_CI_f <- function(fits){
  ## Description
  # Given a list of imputed datasets under the same imputation model
  # it returns the pooled CIs of the regression coefs
  
  ## Example internals
  # fits <- lm_fits[[1]]
  # fits <- m2_mi$blasso
  
  ## Number of multiple imputation
  m <- length(fits)
  
  ## Estiamtes
  coefs <- t(sapply(X = fits,
                    FUN = function(x) coef(x)))
  Q_bar <- colMeans(coefs)
  
  ## Variances
  summa_models <- lapply(X = fits,
                         FUN = function(x) summary(x))
  all_vcov <- lapply(X = summa_models,
                     FUN = function(x) vcov(x))
  U_bar <- diag(Reduce('+', all_vcov) / m)
  
  B <- diag(1 / (m-1) * (t(coefs) - Q_bar) %*% t(t(coefs) - Q_bar))
  
  T <- U_bar + B + B/m
  
  ## Degrees of freedom
  nu <- .miDf(length(fits), b = B, t = T, summa_models[[1]]$df[2])
  
  t_nu <- qt(1 - (1-parms$alphaCI)/2, 
             df = nu)
  
  ## CI computation
  CI <- c(lwr = Q_bar - t_nu * sqrt(T), 
          upr = Q_bar + t_nu * sqrt(T))
  
  return(CI = CI)
}

onetree <- function(xobs, xmis, yobs, s) {
  ## Ripped off mice package mice.impute.rf. Fits one tree for the 
  ## random forest. Used in your own random forest version
  ## based on Doove et al 2014 Algorithm A.1
  
  ## For internals:
  # yobs = y_obs
  # xobs = X_obs
  # xmis = X_mis
  # s = 1
  
  ## Body
  fit <- randomForest::randomForest(x = xobs, y = yobs,
                                    ntree = 1)
    # Takes 1 bootsrap sample, fits 1 tree with random feature selection,
    # provides 1 set of donors.
    # 1 of k repetitions of steps 2a and 2b in algorithm A.1 Doove et al 2014.
    # Each single tree takes 1 bootstrap sample and then goes through
    # partitioning with random sampling of the features. Using this
    # "onetree" function 10 times, you end up with 10 different bootsrapped
    # datasets on which a single tree is fit. The pooling of donors 
    # identified is what makes this a random forest in the end.
  
  # Leaf position of observed part of y (yobs)
  leafnr <- predict(object = fit, newdata = xobs, nodes = TRUE)
  
  # Leaf position of missing part of y (ymis)
  nodes <- predict(object = fit, newdata = xmis, nodes = TRUE)
  
  # Define donors (pool of yobs for a ymis that falls in 1 leaf)
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

find.collinear <- function(x, threshold = 0.999) {
  # Find collinear predictors to be excluded from imputation model
  # Credits to mice package
  nvar <- ncol(x)
  x <- data.matrix(x)
  r <- !is.na(x)
  nr <- apply(r, 2, sum, na.rm = TRUE)
  ord <- order(nr, decreasing = TRUE)
  xo <- x[, ord, drop = FALSE]
  varnames <- dimnames(xo)[[2]]
  z <- suppressWarnings(cor(xo, use = "pairwise.complete.obs"))
  hit <- outer(seq_len(nvar), seq_len(nvar), "<") & (abs(z) >= threshold)
  out <- apply(hit, 2, any, na.rm = TRUE)
  return(varnames[out])
}

# Preprocess Single Imputation
prep_SI <- function(dt_in, 
                    model.var,
                    ...){
  ## Inputs
  # dt_in = Xy_mis
  # iters = 10
  # model.var <- parms$z_m_id
  
  ## Body
  pMat     <- quickpred(dt_in, mincor = .3)
  mids_obj <- mice(dt_in,
                   predictorMatrix = pMat,
                   printFlag       = TRUE,
                   method          = "pmm",
                   ...
                   # ridge           = cond$ridge,
                   # m               = 1,
                   # maxit           = 10,
                   )
  
  dt_out <- complete(mids_obj)
  
  # Define Single Imputed data as the auxiliary set
  target <- which(colnames(dt_in) %in% model.var)
  dt_out[, target] <- dt_in[, target]
  
  # Store outputs
  output <- list(dt_out = dt_out,
                 mids_obj = mids_obj)
  # Return output
  return(output)
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
 
# Rhat convergence checks 
Rhat.sim <- function(out, cond = 1, meth, dat, iter_max = NULL, iter_burn){
  ## Description
  # Given the simulation output with multiple chains, a condition number, 
  # an imputation method name, and a data repetition indicator, it gives 
  # an Rhat qunatifying the ratio of the between and within chain variance
  # The iters_range is used if you want to check convergence at a different
  # point of the chain (other than the complete chain)
  
  ## Internals:
  # out  = out_cnv
  # cond = 1  # some condition
  # meth = out_cnv$parms$method[[1]]
  # dat  = 1
  # iter_max = 250
  # iter_burn = 50
  
  ## Body
  # Check if iter_max is given
  if(is.null(iter_max)) iter_max <- out$parms$iters
    
  # Get the Rhat
  Rhat.sim.out <- sapply(1:out$parms$zm_n, function(v){
    # Extract mean imptued value per variable per iteration
    sims <- sapply(1:out$parms$chains, function(j) {
      rowMeans(out[[dat]][[cond]]$imp_values[[meth]][[j]][[v]][(iter_burn+1):iter_max, ])
    })
    
    # Split chains in half
    half_indx <- 1:floor(nrow(sims)/2)
    split_sims <- cbind(sims[half_indx, ],
                        sims[-half_indx, ])
    
    # Compute Split-Rhat
    Rhat(split_sims)
  })
  names(Rhat.sim.out) <- out$parms$z_m_id
  return(Rhat.sim.out)
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

# Euclidean Distance overall
res_ed_overall <- function(out, condition){
# condition = 1
# out = output
  ## Prep ##
  
  store_avg <- list()
  store_CIC <- list()
  
  for (mod in 1:2) {
    # mod <- 1
    model <- c("m1", "m2")[mod]
    est <- paste0(model, "_EST")
    ci  <- paste0(model, "_CI")
    select_cond <- names(out[[1]])[condition]
    
    ## Step 2. Bias ##
    avg <- sapply(out$parms$methods, function(m){
      # m <- out$parms$methods[2]
      store <- NULL
      count <- 0
      for (i in 1:out$parms$dt_rep) {
        # i <- 3
        succ_method <- colnames(out[[i]][[select_cond]][[est]])
        result <- as.matrix(out[[i]][[select_cond]][[est]])[, succ_method %in% m]
        if(any(is.na(result))) {
          count <- count+1
          next
        }
        store <- cbind(store, result)
        # rownames(store) <- rownames(out[[i]][[select_cond]][[est]])
      }
      c(rowMeans(store, na.rm = TRUE), rep = ncol(store)) # MCMC statistics 
      # rowSums(!is.na(store))
    })
    
    # Store Objects
    avg <- data.frame(avg[, colSums(is.nan(avg)) == 0])
    # get rid of NaNs that come up in exp3 for conditions that are not using 
    # certain methods
    validReps  <- avg["rep", ] # number of successes
    avg <- avg[-which(rownames(avg) == "rep"), ]
    store_avg[[mod]] <- avg
    # store_psd[[mod]] <- avg[, "GS"] # pseudo true values
    
    # storing threshold
    str_thrs <- nrow(out[[1]][[select_cond]][[ci]])/2
    
    # Confidence Interval Coverage
    CIC <- sapply(out$parms$methods, function(m){
      # m <- out$parms$methods[1]
      store <- NULL
      for (i in 1:out$parms$dt_rep) {
        succ_method <- colnames(out[[i]][[select_cond]][[est]])
        col_indx <- succ_method %in% m
        cond_est <- out[[i]][[select_cond]][[est]]
        cond_CI  <- out[[i]][[select_cond]][[ci]]
        ci_low   <- cond_CI[1:str_thrs, ]
        ci_hig   <- cond_CI[-(1:str_thrs), ]
        
        store <- cbind(store, 
                       ci_low[, col_indx] < avg[, "GS"] &
                         avg[, "GS"] < ci_hig[, col_indx]
        )
      }
      rowMeans(store, na.rm = TRUE) # MCMC statistics 
    })
    
    if(is.vector(CIC)){ # if it's a vector make it a dataframe with special care
      CIC_nan <- data.frame(t(is.nan(CIC)))
      CIC <- data.frame(t(CIC[colSums(CIC_nan) == 0]))
    } else {
      CIC <- CIC[, colSums(is.nan(CIC)) == 0] 
      CIC <- data.frame(CIC)
    }
    store_CIC[[mod]] <- CIC
  }

  parms_allMCMC <- do.call(rbind, store_avg)
  parms_psdTrue <- parms_allMCMC[, "GS"]
  
  parms_allCIC <- do.call(rbind, store_CIC)
  
  # Bias
  ed_est <- sapply(parms_allMCMC, function(avg_col){
    dist(rbind(avg_col, parms_psdTrue), method = "euclidean")
  } 
  )
  
  # Confidence intervals
  ed_ci <- sapply(parms_allCIC, function(CIC_col){
    # avg_col <- avg[[1]]
    dist(rbind(CIC_col, rep(.95, nrow(parms_allCIC))), method = "euclidean")
  } 
  )
  return(list(ed_est = data.frame(t(ed_est)),
              ed_ci = data.frame(t(ed_ci)))
  )
}

# Where should this function be loaceted?
res_sum <- function(out, model, condition = 1, bias_sd = FALSE){
  # model = "sem"
  # model = "semR"
  # model = "CFA"
  # model = "m1"
  # condition = 2
  # out = output
  # bias_sd = TRUE
  
  ## Prep ##
  est <- paste0(model, "_EST")
  ci  <- paste0(model, "_CI")
  select_cond <- names(out[[1]])[condition]

  ## Step 2. Bias ##
  # avg <- sapply(out$parms$methods, function(m){
  #   # m <- out$parms$methods[2]
  #   store <- NULL
  #   count <- 0
  #   for (i in 1:out$parms$dt_rep) {
  #     # i <- 1
  #     succ_method <- colnames(out[[i]][[select_cond]][[est]])
  #     result <- as.matrix(out[[i]][[select_cond]][[est]])[, succ_method %in% m]
  #     if(any(is.na(result))) {
  #       count <- count+1
  #       next
  #     }
  #     store <- cbind(store, result)
  #   }
  #   # apply(store, 1, var)
  #   c(rowMeans(store, na.rm = TRUE), rep = ncol(store)) # MCMC statistics 
  # })
  
  par_avg <- NULL
  par_var <- NULL
  for(m in out$parms$methods){
    # m <- out$parms$methods[2]
    store <- NULL
    count <- 0
    for (i in 1:out$parms$dt_rep) {
      # i <- 1
      succ_method <- colnames(out[[i]][[select_cond]][[est]])
      result <- as.matrix(out[[i]][[select_cond]][[est]])[, succ_method %in% m]
      if(any(is.na(result))) {
        count <- count+1
        next
      }
      store <- cbind(store, result)
    }
    par_avg <- cbind(par_avg, 
                     c(rowMeans(store, na.rm = TRUE), 
                       rep = ncol(store)))
    par_var <- cbind(par_var, 
                     apply(store, 1, var))
  }
  colnames(par_avg) <- colnames(par_var) <- out$parms$methods
  
  # Store Objects
  avg <- par_avg
  avg <- data.frame(avg[, colSums(is.nan(avg)) == 0])
    # get rid of NaNs that come up in exp3 for conditions that are not using 
    # certain methods
  validReps  <- avg["rep", ] # number of successes
  avg        <- avg[-which(rownames(avg) == "rep"), ]
  psd_tr_vec <- avg[, "GS"] # pseudo true values
  
  if(out$parms$exp == 1 & model == "sem"){
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
                    bias / psd_tr_vec * 100)
  
  meths <- out$parms$methods[-which(out$parms$methods == "GS")]
  # Bias Mean Standardized
  if(bias_sd == TRUE){
    # Empirical Standard error
    sd_emp <- sapply(out$parms$methods, function(m){
      # m <- out$parms$methods[1]
      store <- NULL
      for (i in 1:out$parms$dt_rep) {
        succ_method <- colnames(out[[i]][[select_cond]][[est]])
        store <- cbind(store, 
                       out[[i]][[select_cond]][[est]][,
                                                      succ_method %in% m]
        )
      }
      return( apply(t(store), 2, sd, na.rm = TRUE) )
    })
    
    bias_sd_out <- (avg - psd_tr_vec) / sd_emp
    
  } else {
    bias_sd_out <- NULL
  }
  
  ## Step 3. CI Coverage ##
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
      ci_low   <- cond_CI[1:str_thrs, col_indx]
      ci_hig   <- cond_CI[-(1:str_thrs), col_indx]
      store <- cbind(store, 
                     ci_low < psd_tr_vec & psd_tr_vec < ci_hig
      )
    }
    rowMeans(store, na.rm = TRUE) # MCMC statistics 
  })
  
  if(is.vector(CIC)){ # if it's a vector make it a dataframe with special care
    CIC_nan <- data.frame(t(is.nan(CIC)))
    CIC <- data.frame(t(CIC[colSums(CIC_nan) == 0]))
  } else {
    CIC <- CIC[, colSums(is.nan(CIC)) == 0] 
    CIC <- data.frame(CIC)
  }
    # get rid of NaNs that come up in exp3 for conditions that are not using 
    # certain methods
  rownames(CIC) <- rownames(bias)
  
  ## Step 4 - Euclidean distances
  # MCMC estimates (per model)
  ed_est <- sapply(avg, function(avg_col){
    # avg_col <- avg[[1]]
    dist(rbind(avg_col, psd_tr_vec), method = "euclidean")
  } 
  )
  
  # Confidence intervals
  ed_ci <- sapply(CIC, function(CIC_col){
    # avg_col <- avg[[1]]
    dist(rbind(CIC_col, rep(.95, nrow(CIC))), method = "euclidean")
  } 
  )
  
  ## Step 5 - Confidence interval width
  # Confidence Interval Coverage
  CIW <- sapply(out$parms$methods, function(m){
    store <- NULL
    for (i in 1:out$parms$dt_rep) {
      succ_method <- colnames(out[[i]][[select_cond]][[est]])
      col_indx <- succ_method %in% m
      cond_est <- out[[i]][[select_cond]][[est]]
      cond_CI  <- out[[i]][[select_cond]][[ci]]
      ci_low   <- cond_CI[1:str_thrs, ]
      ci_hig   <- cond_CI[-(1:str_thrs), ]
      
      store <- cbind(store, (ci_hig[, col_indx] - ci_low[, col_indx]))
    }
    rowMeans(store, na.rm = TRUE) # MCMC statistics 
  })
  CIW <- CIW[, colSums(is.nan(CIW)) == 0] 
  CIW <- data.frame(CIW)
  rownames(CIW) <- rownames(bias)
  
  # Output
  results <- list(cond      = select_cond,
                  MCMC_est  = round(cbind(ref=psd_tr_vec, avg), 3),
                  bias_raw  = round(cbind(ref=psd_tr_vec, bias), 3),
                  bias_per  = bias_per,
                  bias_sd   = bias_sd_out,
                  ci_cov    = round(CIC*100, 3),
                  ed_est    = data.frame(t(ed_est)),
                  ed_ci     = data.frame(t(ed_ci)),
                  CIW       = CIW,
                  var_est   = par_var,
                  validReps = validReps)
  return(results)
}

# Convergence checks
mean_traceplot <- function(out, 
                           dat = 1, # which data repetition
                           method = "blasso", # same name as in parms
                           y_range = c(-10, 20),
                           v_range = NULL,
                           iters = 1:5){
  ## Internals
  # dat = 1
  # out <- out_DURR
  # iters = 1:7
  # v_range = 170:180
  # method = "DURR_all" # same name as in parms
  
  ## Description
  # It prints the traceplots for the mean imputed values in each iteration
  # in different chains, by variable, one dataset, one imputation method
  
  # Variables to display
  if(is.null(v_range)){
    v_range <- 1:length(out$parms$z_m_id)
  }
  
  # Display in same pane
  par(mfrow = c(3, ceiling(length(v_range)/3)))
  # Plot
  # Are imputations of mids class?
  if(class(out[[dat]][[1]]$imp_values[[method]]) == "mids"){
    plot(out[[dat]][[1]]$imp_values[[method]],
         y = paste0("z", v_range))
  } else {
    # For each variable imputed
    for (v in v_range) {
      # v <- v_range[1]
      # CHAIN 1
      # Mean imputed value (across individuals), in the first CHAIN,
      # at each iteration
      mean_imp <- rowMeans(out[[dat]][[1]]$imp_values[[method]][[1]][[v]][iters, ])
      
      plot(iters, mean_imp, type = "l",
           main = method,
           ylim = y_range,
           ylab = paste0("z", v), xlab = "Iteration")
      
      # CHAIN 2 to m 
      # Mean imputed value (across individuals), in the first CHAIN,
      # at each iteration
      for (i in 2:(out$parms$chains)) {
        mean_imp <- rowMeans(out[[dat]][[1]]$imp_values[[method]][[i]][[v]][iters, ])
        lines(iters, mean_imp)
      }
    }
  }
}

# Crossvalidation
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
    # Within the same condition, take mean of average fmi from 
    # each data repetition
    store_0[[i]] <- rowMeans(store_1)
  }
  ridge_range    <- length(unique(out$conds$ridge))
  names(store_0) <- rep(unique(out$conds$ridge), nrow(out$conds)/ridge_range)
  
  # Within the same condition, take mean of the m1 and m2 fmis
  # (previously aggregated across datasets)
  avg_fmi <- round(sapply(store_0, mean), 3)
  
  # Select ridge value with lowest average FMI
  i <- 1; j <- i + ridge_range - 1
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

# Plot function for experiment 1 and 2
plot_fg <- function(dt,
                    type = "bias",
                    parPlot,
                    dt_reps = 500,
                    ci_lvl = .95,
                    axis.name.x = NULL,
                    cond_labels = NULL,
                    plot_name = NULL,
                    y_axLab = TRUE,
                    summy = FALSE,
                    meth_compare) {
  ## Function inputs
  ## Generic
  # type = "bias"
  # axis.name.x = NULL #"Axis Name Here"
  # plot_name = NULL #"Untitled"
  # ci_lvl = .95 # confidence interval level
  # dt_reps = 500 # MCMC repetitions
  # y_axLab = TRUE # say I want the labels
  # parPlot <- list(means = 1:6,
  #                 variances = 7:12,
  #                 covariances = 13:27)
  # summy = TRUE # requires summary or not of stats
  # meth_compare = rev(c("DURR_la", "IURR_la", "blasso",# "bridge",
  #                      "MI_PCA",
  #                      "MI_CART", "MI_RF", "missFor", "CC"))
  # cond_labels = NULL
  ## EXP 1
  # dt = lapply(1:length(res$sem),
  #             function(x) data.frame( res$sem[[x]]$bias_per))
  # dt = lapply(1:length(res$sem),
  #             function(x) data.frame( res$sem[[x]]$ci_cov))
  ## EXP 2
  # dt = lapply(1:length(res$semR),
  #             function(x) data.frame( res$semR[[x]]$bias_per))[1:4]
  # dt = lapply(1:length(res$semR),
  #             function(x) data.frame( res$semR[[x]]$bias_sd))[1:4]
  # dt = lapply(1:length(res$semR),
  #             function(x) data.frame( res$semR[[x]]$ci_cov))[1:4]
  # dt = lapply(1:length(res$CFA),
  #             function(x) data.frame( res$CFA[[x]]$bias_per))[1:4]
  # parPlot = list(Loadings = 1:10)
  # means = 1:10
  # varis = 11:16
  # covas = 21:65
  ## EXP 5
  # dt = lapply(1:length(res$semR),
  #             function(x) data.frame( res$semR[[x]]$bias_per))[1:4]
  # meth_compare = c("DURR_all","DURR_si","IURR_all","IURR_si","blasso",
  #                  #"bridge",
  #                  "MI_PCA","MI_CART" ,"MI_RF","MI_OP",
  #                  "missFor","CC")
  ## Prep data for plot
  dt_preEdit <- lapply(parPlot, function(x){
    lapply(dt, function(d){
      d[x, meth_compare]
    })
  })
  
  # Apply if summary version required
  if(summy == TRUE){
  dt_preEdit <- lapply(parPlot, function(x){
    # x <- parPlot[[3]]
    lapply(dt, function(d){
      # d <- dt[[2]]
      temp <- abs(d[x, meth_compare])
      as.data.frame(
        sapply(temp, function(j){ summary(j)[c(1, 4, 6)] })
      )
    })
  })
  }
  
  n_parms <- length(parPlot) # number of parameters
  n_conds <- length(dt_preEdit[[1]]) # number of conditions
  n_situa <- n_conds * n_parms # number of combinations/situations/plots
  
  # Bring list to 1 level
  dt_edit <- Reduce(c, dt_preEdit)

  # Count the number of estimates (for each situation)
  n_estXcond <- lapply(dt_edit, nrow)
  
  # Count the number of estimates (for each parameter type)
  n_estXparm <- unlist(n_estXcond[seq(1, n_situa, n_conds)])

  # Make names prettier
  dt_edit <- lapply(dt_edit, function(x){
    colnames(x) <- sub("_la", "", colnames(x))
    colnames(x) <- sub("_", "-", colnames(x))
    return(x)
  })
  
  # Add Blank Row to improve readability
  dt_edit <- lapply(dt_edit, function(x){
    x[nrow(x)+1, ] <- 0  # add blank row to improve visualization  
    return(x)
  })
  
  # Count contents
  n <- lapply(dt_edit, nrow)
  
  # Shape for ggplot
  dt_edit <- lapply(dt_edit, gather)
  dt_edit <- lapply(dt_edit, function(x){
    x$id <- 1:nrow(x)
    return(x)
  })
  n_facet <- length(dt_edit)
  
  # Combine for facet
  dt_edit <- do.call(rbind, dt_edit)
  
  plot_step <- length(dt)
  gpf1 <- seq(0, length(dt)*length(parPlot), by = plot_step)
  
  # Grid Plot Factor 1
  conds_list <- lapply(gpf1[-1], function(x){
    rep(1:length(dt), each = n[[x]]*length(meth_compare))
  })
  conds <- do.call(c, conds_list)
  
  # Grid Plot Factor 2
  # parT_index <- (sapply(parPlot, length)+1)*length(dt)*length(meth_compare)
  parT_index <- (n_estXparm + 1) * length(dt)*length(meth_compare)
  parT <- lapply(1:length(parPlot), function(x){
    rep(x, parT_index[x])
  })
  parT <- Reduce(c, parT)
  
  # Final Data prep
  if(is.null(cond_labels)){cond_labels <- 1:length(dt)}
  dt_edit <- cbind(dt_edit, 
                   conds = factor(cond_labels[conds], levels = cond_labels),
                   parT = factor(parT, labels = names(parPlot)))
  
  # Define Step Size for all parameters sets
  step_size <- vector("list", 3)
  plot_steps <- vector("list", 3)
  plot_ybreaks <- vector("list", 3)
  plot_hlines <- vector("list", 3)
  for(i in 1:n_parms){
    ref_data <- dt_preEdit[[i]][[1]]
    step_size[[i]]   <- (
      nrow(ref_data) + # number of parameters
      1 # account for additional empty row
    ) / 2 # place label in the middle
    plot_steps[[i]] <- seq(0, 
                        (nrow(ref_data)+1) *  # number of parameters + empty row
                        ncol(ref_data), # number of methods
                        by = step_size[[i]])
    plot_ybreaks[[i]] <- plot_steps[[i]][c(FALSE, TRUE)]    # position label skip
    plot_hlines[[i]] <- plot_steps[[i]][c(TRUE, FALSE)]  
  }
  
  # Methods labels
  plot_ylabels <- as.character(unique(dt_edit$key)) # unique for everyone
  
  # Ticks 
  if(type == "bias"){
    # Plot Limits
    plot_xlim   <- c(-20, 20)
    if(summy == TRUE){
      plot_xlim <- c(0, 40)
    }
    # X axis
    xbreaks_low <- min(plot_xlim)
    xbreaks_top <- max(plot_xlim)
    xbreaks_center <- mean(plot_xlim)
    xbreaks_midlow <- mean(xbreaks_center:xbreaks_low)
    xbreaks_midtop <- mean(xbreaks_center:xbreaks_top)
    plot_xbreaks <- c(xbreaks_low, xbreaks_midlow, 
                      xbreaks_center, 
                      xbreaks_midtop, xbreaks_top)# c(-20, -10, 0, 10, 20)
    plot_xlabels <- as.character(plot_xbreaks)#c("-20", "-10", "0", "10", "20")
    plot_vlines <- c(xbreaks_midlow, xbreaks_midtop)
    if(summy == TRUE){
      plot_vlines <- c(xbreaks_midlow)
    }
  }
  if(type == "ci"){
    # Redefine values as differences from target
    dt_edit$value[dt_edit$value != 0] <- dt_edit$value[dt_edit$value != 0] - 95
    
    # Plot Limits (reference: 0 = .95)
    plot_xlim   <- c(-10, 5)
    
    # SE for threshold 
    SEp <- sqrt(ci_lvl*(1-ci_lvl)/dt_reps)
    low_thr <- ((.95-SEp*2)-.95)*100
    hig_thr <- ((.95+SEp*2)-.95)*100
    
    # X axis
    plot_xbreaks <- c(min(plot_xlim), -5, low_thr, 0, hig_thr, max(plot_xlim))
    plot_xlabels <- as.character(round((plot_xbreaks+95)/100, 2))
    plot_vlines <- c(-5, low_thr, hig_thr, 4.99)
  }
  
  # Colors and texts
  font.plot        <- "Arial" # font for the whole plot
  x.axis.text.size <- 7.5 # Scale of plotted numbers
  y.axis.text.size <- 7.5 # Imputation Methods names
  grid.text.size   <- 7.5 # Condition + Parameter type
  segme.thick      <- .5 # thickness of lines reporting results
  segme.color      <- "black" # color  of lines reporting results
  h.lines.thick    <- .10 # thickness of lines separating methods
  h.lines.color    <- "black" # color of lines separating methods
  v.lines.thick    <- .375 # thickness of reference lines
  v.lines.color    <- "darkgray" # color of reference lines
  v.lines.type     <- "dashed" # line type of reference lines
  
  # Plot
  p <- ggplot(dt_edit, aes(x = value, y = id)) +
    # Title and axis labels
    labs(title = plot_name,
         x     = axis.name.x, 
         y     = element_blank()) +
    theme(plot.title   = element_text(#family = font.plot,
                                      size = 6.5, 
                                      face = "plain", 
                                      hjust = .5,
                                      vjust = .5),
          axis.title.x = element_text(size = 5,
                                      face = "plain",
                                      family = font.plot),
          # Scale of plotted numbers
          axis.text.x  = element_text(#family = font.plot,
                                      size = 5,
                                      angle = 90,
                                      vjust = .5),
          # Imputation Methods names
          axis.text.y  = element_text(#family = font.plot,
                                      size = y.axis.text.size),
          plot.margin  = unit(c(0, .0, .0, .0), "cm"),
          # Background
          panel.background = element_rect(fill = "white", colour = "white"),
          panel.border     = element_rect(colour = "lightgray", fill = NA),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          # Condition + Parameter type
          strip.text = element_text(#family = font.plot,
                                    size = grid.text.size,
                                    face = "plain",
                                    margin = unit(c(.10, .10, .10, .10), "cm")) 
    )
  
  # Plot Content
  if(length(levels(dt_edit$parT)) == 1){
    # For Factor Loadings
    p <- p +
      geom_segment(aes(xend = 0, 
                       yend = id),
                   data = dt_edit,
                   size = segme.thick,
                   color = segme.color)
  }
  if(length(levels(dt_edit$parT)) > 2){
    p <- p +
    # For means and variances
      geom_segment(aes(xend = 0, 
                       yend = id),
                   data = dt_edit[dt_edit$parT != "covariances", ],
                   size = segme.thick,
                   color = segme.color) +
    # For covariances (lighter tone)
      geom_segment(aes(xend = 0, 
                       yend = id),
                   data = dt_edit[dt_edit$parT == "covariances", ],
                   size = segme.thick,
                   color = segme.color)  
  }
  # Axis
  p <- p +
    # X Axis
    scale_x_continuous(breaks = plot_xbreaks,
                       labels = plot_xlabels) +
    geom_vline(xintercept = plot_vlines,
               size       = v.lines.thick,
               linetype   = v.lines.type, 
               color      = v.lines.color) +
    coord_cartesian(xlim  = plot_xlim)
  
  # Add Facet and y lines
  if(length(parPlot) == 1){
    p <- p + 
      facet_grid_custom(rows = vars(parT),
                        cols = vars(conds),
                        scales = "free", 
                        scale_overrides = list(
                          scale_override(1, scale_y_continuous(breaks = plot_ybreaks[[1]],
                                                               labels = plot_ylabels))
                        )) +   
      # Horizontal lines (method separation)
      geom_hline(data = data.frame(yint = plot_hlines[[1]], 
                                   parT = levels(dt_edit$parT)[[1]]),
                 aes(yintercept = yint), 
                 size = h.lines.thick, 
                 color = h.lines.color)
  }
  if(length(parPlot) == 2){
    p <- p + 
      facet_grid_custom(rows = vars(parT),
                        cols = vars(conds),
                        scales = "free", 
                        scale_overrides = list(
                          scale_override(1, scale_y_continuous(breaks = plot_ybreaks[[1]],
                                                               labels = plot_ylabels)),
                          scale_override(2, scale_y_continuous(breaks = plot_ybreaks[[2]],
                                                               labels = plot_ylabels))
                        )) +   
      # Horizontal lines (method separation)
      geom_hline(data = data.frame(yint = plot_hlines[[1]], 
                                   parT = levels(dt_edit$parT)[[1]]),
                 aes(yintercept = yint), 
                 size = h.lines.thick, 
                 color = h.lines.color) +
      geom_hline(data = data.frame(yint = plot_hlines[[2]], 
                                   parT = levels(dt_edit$parT)[[2]]),
                 aes(yintercept = yint), 
                 size = h.lines.thick, 
                 color = h.lines.color)
  }
  if(length(parPlot) == 3){
    p <- p + 
      facet_grid_custom(rows = vars(parT),
                        cols = vars(conds),
                        scales = "free", 
                        scale_overrides = list(
                          scale_override(1, scale_y_continuous(breaks = plot_ybreaks[[1]],
                                                               labels = plot_ylabels)),
                          scale_override(2, scale_y_continuous(breaks = plot_ybreaks[[2]],
                                                               labels = plot_ylabels)),
                          scale_override(3, scale_y_continuous(breaks = plot_ybreaks[[3]],
                                                               labels = plot_ylabels))
                        )) +   
      # Horizontal lines (method separation)
      geom_hline(data = data.frame(yint = plot_hlines[[1]], 
                                   parT = levels(dt_edit$parT)[[1]]),
                 aes(yintercept = yint), 
                 size = h.lines.thick, 
                 color = h.lines.color) +
      geom_hline(data = data.frame(yint = plot_hlines[[2]], 
                                   parT = levels(dt_edit$parT)[[2]]),
                 aes(yintercept = yint), 
                 size = h.lines.thick, 
                 color = h.lines.color) +
      geom_hline(data = data.frame(yint = plot_hlines[[3]], 
                                   parT = levels(dt_edit$parT)[[3]]),
                 aes(yintercept = yint), 
                 size = h.lines.thick, 
                 color = h.lines.color)
  }
  
  # Visualize Plot
  p
  return(p)
}

# Custom plot_gg function for experiment 4
plot_exp4 <- function(dt, 
                      dt_CIW = NULL,
                      dt_reps = 500,
                      ci_lvl = .95,
                      type = "ci", 
                      plot_cond = NULL,
                      plot_name = NULL,
                      bar_col = "#595959",
                      meth_compare,
                      meth_sort = FALSE) {
  
  # Function Inputs
  ## New Input
  # dt = list(lapply(1:length(res$m1),
  #                  function(x) res$m1[[x]]$bias_per["rel", ]),
  #           lapply(1:length(res$m2),
  #                  function(x) res$m2[[x]]$bias_per["NatAt", ]))
  # # CI Par interest
  # dt = list(lapply(1:length(res$m1),
  #                  function(x) res$m1[[x]]$ci_cov["rel", ]),
  #           lapply(1:length(res$m2),
  #                  function(x) res$m2[[x]]$ci_cov["NatAt", ]))
  # # Mean CIW
  # dt = list(lapply(1:length(res$m1),
  #                  function(x) data.frame(t(colMeans(res$m1[[x]]$CIW)))),
  #           lapply(1:length(res$m2),
  #                  function(x) data.frame(t(colMeans(res$m2[[x]]$CIW)))))
  # # ED
  # dt = list(lapply(1:length(res$m1),
  #                  function(x) res$m1[[x]]$ed_est),
  #           lapply(1:length(res$m2),
  #                  function(x) res$m2[[x]]$ed_est))
  # # CIED
  # dt = list(lapply(1:length(res$m1),
  #                  function(x) res$m1[[x]]$ed_ci),
  #           lapply(1:length(res$m2),
  #                  function(x) res$m2[[x]]$ed_ci))
  ## Generic inputs
  # dt_reps = 500
  # ci_lvl = .95
  # type = c("bias", "ci", "ciw", "ed")[1]
  # dt_CIW = NULL
  # plot_name = "Untitled"
  # plot_cond = "(empty)"
  # meth_compare = c("DURR_la", "IURR_la", "blasso", "bridge",
  #                  "MI_PCA",
  #                  "MI_CART", "MI_RF", "missFor", "CC")
  # meth_sort = FALSE
  # bar_col = "darkgray"
  
  # Gather data within list
  dt_edit <- lapply(Reduce(c, dt), gather)
  
  # Clean methods
  dt_edit <- lapply(dt_edit, function(x){
    # Make methods rownames for ordering reasons
    rownames(x) <- x$key
    # select methods to keep
    x <- x[meth_compare, ] 
    # Make names prettier
    x$key <- sub("_la", "", x$key)
    x$key <- sub("_", "-", x$key)
    # fix order of methods
    x$key <- factor(x$key, levels = x$key)
    return(x)
  })
  
  # Add Condition tag
  dt_edit <- lapply(1:4, function(x){
    if(x<3){
      dt_edit[[x]]$cond <- paste0("Condition ", x)
    } else {
      dt_edit[[x]]$cond <- paste0("Condition ", x-2)
    }
    return(dt_edit[[x]])
  })
  
  # Add model tag
  dt_edit <- lapply(1:4, function(x){
    if(x<3){
      dt_edit[[x]]$model <- "Model 1"
    } else {
      dt_edit[[x]]$model <- "Model 2"
    }
    return(dt_edit[[x]])
  })
  
  # Gather data for plot
  dt_edit <- Reduce(rbind, dt_edit)
  
  # Summary specific
  if(type == "bias"){
    plot_xlim   <- c(-.20, .20) * 100
    plot_breaks <- c(-.20, -.1, 0, .1, .20)  * 100
    plot_labels <- c("-20", "-10", "0", "10", "20")
    plot_vlines <- c(-.1, .1)  * 100
    plot_hlines <- 1:length(meth_compare)
  }
  
  if(type == "ci"){
    dt_edit$value <- dt_edit$value - 95    # Transform to a difference value
    plot_xlim   <- c(-5, 5)                # Plot Limits (reference: 0 = .95)
    SEp <- sqrt(ci_lvl*(1-ci_lvl)/dt_reps) # SE for threshold 
    low_thr <- ((.95-SEp*2)-.95)*100
    hig_thr <- ((.95+SEp*2)-.95)*100
    plot_breaks <-  c(-5, low_thr, 0, hig_thr, 5)
    plot_labels <- as.character(round((plot_breaks+95)/100, 2))
    plot_hlines <- 1:length(meth_compare)
    plot_vlines <- c(-5, low_thr, hig_thr)
  }
  
  if(type == "ed"){
    plot_xlim   <- c(0, round_any(max(dt_edit[dt_edit$key != "bridge", "value"]),  
                                  0.1, 
                                  f = ceiling)
    )
    plot_breaks <- round_any(seq(plot_xlim[1], plot_xlim[2], 
                                 length.out = 3), 0.1, f = ceiling)
    plot_breaks <- round(seq(plot_xlim[1], plot_xlim[2], length.out = 3), 
                         1)
    plot_labels <- gsub("0.", ".", plot_breaks)
    # if("GS" %in% dt_edit$key) {
    #   plot_vlines <- dt_edit[dt_edit$key=="GS", "value"]
    # } else {
      plot_vlines <- NULL
    # }
    plot_hlines <- NULL
  }
  
  if(type == "ciw"){
    ref_max <- max(dt_edit[dt_edit$key == "CC", "value"])
    ref_min <- min(dt_edit[dt_edit$key == "missFor", "value"])
    plot_xlim   <- c(0, (ref_max*1.5))
    plot_breaks <- c(0, ref_min, ref_max, ref_max*1.5)
    plot_labels <- as.character(round(plot_breaks, 2))
    plot_hlines <- 1:length(meth_compare)
    plot_vlines <- NULL
  }
  
  # Plot
  p <- ggplot(dt_edit, aes(x = value, y = key)) +
    # Content
    geom_bar(stat = "identity",
             fill = bar_col,
             position = "dodge",
             width = .5) +
    # Faceting
    facet_grid(rows = vars(model),
               cols = vars(cond),
               scales = "free") +
    # Cosmetic
    labs(title = element_blank(),
         x     = element_blank(), 
         y     = element_blank()) + 
    theme(axis.text.x  = element_text(size = 5,
                                      angle = 90,
                                      vjust = .5),
          axis.text.y  = element_text(size = 5),
          plot.margin  = unit(c(.05, .0, .0, .0), "cm"),
          # Background
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          # Facet Related
          strip.text = element_text(size = 10,
                                    face = "plain",
                                    margin = unit(c(.10, .10, .10, .10), "cm")) 
    ) + 
    # X Axis
    scale_x_continuous(breaks = plot_breaks,
                       labels = plot_labels) +
    coord_cartesian(xlim = plot_xlim) + 
    
    # Horizontal lines
    geom_hline(yintercept = plot_hlines+.5, 
               size = .25,
               color = "gray") +
    
    # Vertical lines
    geom_vline(xintercept = plot_vlines,
               size = .375,
               linetype = "dashed", 
               color = "black")
  p
  return(p)
}

plot_exp4_meth <- function(dt,
                    type = "bias",
                    dt_reps = 500,
                    ci_lvl = .95,
                    focal = "", # name of prameter to highlight
                    small.ef = "",
                    meth_compare) {
  ## Function inputs
  ## Generic
  # type = "bias"
  # ci_lvl = .95 # confidence interval level
  # dt_reps = 500 # MCMC repetitions
  # y_axLab = TRUE # say I want the labels
  # meth_compare = rev(c("DURR_la", "IURR_la", "blasso", "bridge",
  #                      "MI_PCA",
  #                      "MI_CART", "MI_RF", "missFor", "CC", "MI_OP"))
  # focal = "rel"
  # small.ef = ""
  # dt = lapply(1:length(res$m1),
  #             function(x) res$m1[[x]]$bias_per)
  # dt = lapply(1:length(res$m2),
  #             function(x) res$m2[[x]]$bias_per)
  # dt = lapply(1:length(res$m1),
  #             function(x) res$m1[[x]]$bias_raw)
  # dt = lapply(1:length(res$m1),
  #             function(x) res$m1[[x]]$bias_sd)
  # dt = lapply(1:length(res$m1),
  #             function(x) res$m1[[x]]$ci_cov)
  # dt = lapply(1:length(res$m1),
  #             function(x) res$m1[[x]]$CIW)
  
  ## Prep data for plot (take absolute value)
  dt_preEdit <- lapply(dt, function(d){
    abs(d[, meth_compare])
  })
  
  # conds <- rep(c("Condition 1", "Condition 2"), 2)
  conds <- rep(c("low-dim \n p = 243, n = 1000", 
                 "high-dim \n p = 243, n = 300"), 2)
  
  # CONTINUE FROM HERE
  dt_edit <- lapply(1:length(dt_preEdit), function(id) {
    x <- dt_preEdit[[id]]
    
    # Make Methods names prettier
    colnames(x) <- sub("_la", "", colnames(x))
    colnames(x) <- sub("_", "-", colnames(x))
    
    # Extract Methods Name
    methods <- names(x)
    
    # Extract Results for a method
    output_2 <- lapply(1:ncol(x), function(l){
      par_names <- rownames(x[l])[order(x[l], decreasing = TRUE)]
      par_value <- x[l][order(x[l], decreasing = TRUE), ]

      # Compose output for method
      output_1 <- data.frame(key = methods[l], # method
                             par = c(par_names, ""),
                             value = c(par_value, 0),
                             conds = factor(conds[id]))
      return(output_1)
    })
    
    return(output_2)
  })
  
  # Put them in groups by condition
  dt_edit <- lapply(dt_edit, function(x){do.call(rbind, x)})
  
  # Count contents
  n <- lapply(dt_preEdit, function(x){nrow(x)+1}) # +1 to improve spacing
  
  dt_edit <- lapply(dt_edit, function(x){
    x$id <- 1:nrow(x)
    return(x)
  })
  n_facet <- length(dt_edit)
  
  # Combine for facet in one list
  dt_edit <- do.call(rbind, dt_edit)
  
  # Define Step Size for all parameters sets
  step_size   <- (
    nrow(dt_preEdit[[1]]) + # number of parameters
      1 # account for additional empty row
  ) / 2 # place label in the middle
  plot_steps <- seq(0, 
                    (nrow(dt_preEdit[[1]])+1) *  # number of parameters + empty row
                      ncol(dt_preEdit[[1]]), # number of methods
                    by = step_size)
  plot_ybreaks <- plot_steps[c(FALSE, TRUE)]    # position label skip
  plot_hlines <- plot_steps[c(TRUE, FALSE)]  
  
  # Methods labels
  plot_ylabels <- as.character(unique(dt_edit$key)) # unique for everyone
  
  # Paramter Labels
  # dt_edit$par[!dt_edit$par %in% focal] <- ""
  # Ticks 
  if(type == "bias"){
    # Levels order
    levs <- c(yes = ">10%", no = "<10%", focal = "Focal")
    
    # Grid Plot Color based on exceeding or not PRB reference
    flag <- ifelse(dt_edit$value >= 10, yes = levs[1], no = levs[2])
    flag[dt_edit$par %in% focal] <- levs[3]
    # flag[dt_edit$par %in% small.ef] <- "Largest Bias"
    dt_edit$flag <- factor(flag, levels = levs)
    
    # Plot Limits
    plot_xlim   <- c(0, 70)
    
    # X axis
    plot_xbreaks <- seq(min(plot_xlim), max(plot_xlim), by = 10)
    plot_xlabels <- as.character(plot_xbreaks)
    plot_vlines <- c(10)
  }

  if(type == "ci"){
    # Redefine values as differences from target
    dt_edit$value[dt_edit$value != 0] <- dt_edit$value[dt_edit$value != 0] - 95
    
    # Plot Limits (reference: 0 = .95)
    plot_xlim   <- c(-5, 5)
    
    # SE for threshold 
    SEp <- sqrt(ci_lvl*(1-ci_lvl)/dt_reps)
    low_thr <- ((.95-SEp*2)-.95)*100
    hig_thr <- ((.95+SEp*2)-.95)*100
    
    # Levels order
    levs <- c(yes = "Significant", 
              no = "Not significant", 
              focal = "Focal")
    
    # Grid Plot Color based on exceeding or not PRB reference
    flag <- ifelse(dt_edit$value >= hig_thr | dt_edit$value <= low_thr, 
                     yes = levs[1], no = levs[2])
    flag[dt_edit$par %in% focal] <- levs[3]
    dt_edit$flag <- factor(flag, levels = levs)
    
    # X axis
    plot_xbreaks <- c(-5, low_thr, 0, hig_thr, 5)
    plot_xlabels <- as.character(round((plot_xbreaks+95)/100, 2))
    plot_vlines <- c(-5, low_thr, hig_thr)
  }

  # Colors and texts
  font.plot        <- "Arial" # font for the whole plot
  x.axis.text.size <- 7.5 # Scale of plotted numbers
  y.axis.text.size <- 7.5 # Imputation Methods names
  grid.text.size   <- 10 # Condition + Parameter type
  segme.thick      <- 1 # thickness of lines reporting results (was 1)
  small.color      <- "darkgray" # color of lines |PRB| < 10% 
  large.color      <- "lightgray" # color of lines |PRB| < 10% 
  focal.color      <- "black" # color of line reporting focal parameter
  h.lines.thick    <- .10 # thickness of lines separating methods (was .375)
  h.lines.color    <- "black" # color of lines separating methods (was gray)
  v.lines.thick    <- .375 # thickness of reference lines
  v.lines.color    <- "darkgray" # color of reference lines
  v.lines.type     <- "dashed" # line type of reference lines
  
  # Plot
  p <- ggplot(dt_edit, aes(x = value, y = id)) +
    # Title and axis labels
    theme(plot.title   = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x  = element_text(size = x.axis.text.size),
          axis.text.y  = element_text(size = y.axis.text.size),
          plot.margin  = unit(c(.05, .0, .0, .0), "cm"),
          # Background
          panel.background = element_rect(fill = "white", colour = "white"),
          panel.border     = element_rect(colour = "lightgray", fill = NA),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          # Condition + Parameter type (Facet Related)
          strip.text = element_text(#family = font.plot,
            size = grid.text.size,
            face = "plain",
            margin = unit(c(.10, .10, .10, .10), "cm")),
          # Legend
          legend.text = element_text(size = 8),
          legend.position = "bottom"
    ) +
    
    # Plot Content
    geom_segment(aes(xend = 0, 
                     yend = id,
                     colour = flag),
                 size = segme.thick) +
    scale_color_manual(name = "",
                       values = c(large.color,
                                  small.color, 
                                  focal.color)) +
    
    # X Axis
    scale_x_continuous(breaks = plot_xbreaks,
                       labels = plot_xlabels) +
    geom_vline(xintercept = plot_vlines,
               size = v.lines.thick,
               linetype = v.lines.type, 
               color = v.lines.color) +
    coord_cartesian(xlim = plot_xlim) + 
    
    # Y axis
    scale_y_continuous(breaks = plot_ybreaks,
                       labels = plot_ylabels) +
    geom_hline(yintercept = plot_hlines,
               size = h.lines.thick,
               color = h.lines.color) + 
    
    # Facet
    facet_grid(cols = vars(conds))
  
  p
  return(p)
}

plot_time <- function(dt, 
                      plot_cond = NULL,
                      plot_name = NULL,
                      bar_col = "#595959",
                      meth_compare,
                      meth_sort = FALSE) {
  
  # Function Inputs
  ## New Input
  # out_time <- sapply(1:length(names(res[[1]])), res_sem_time, out = res)
  # colnames(out_time) <- names(res[[1]])
  # dt = t(out_time)
  ## Generic inputs
  # dt_reps = 500
  # ci_lvl = .95
  # type = c("bias", "ci", "ciw", "ed")[3]
  # dt_CIW = NULL
  # plot_name = "Untitled"
  # plot_cond = "(empty)"
  # meth_compare = c("DURR_la", "IURR_la", "blasso", "bridge",
  #                  "MI_PCA",
  #                  "MI_CART", "MI_RF", "missFor", "CC")
  # meth_sort = FALSE
  # bar_col = "#595959"
  
  # Gather data within list
  dt_edit <- gather(as.data.frame(dt))
  
  # Clean methods
  # select methods to keep
  dt_edit <- dt_edit[dt_edit$key %in% meth_compare, ] 
  # Make names prettier
  dt_edit$key <- sub("_la", "", dt_edit$key)
  dt_edit$key <- sub("_", "-", dt_edit$key)
  # fix order of methods
  dt_edit$key <- factor(dt_edit$key, levels = rev(unique(dt_edit$key)))
  
  # Add Condition tag
  n_conds <- nrow(dt)
  dt_edit$cond <- paste0("Condition ", 
                         rep(1:n_conds, 
                             nrow(dt_edit)/n_conds))
  # Plot limits
  plot_xlim   <- c(0, 90)
  plot_breaks <- c(0, 30, 60, 90)
  plot_labels <- c("", "30min", "1h", "1h 30min")
  plot_vlines <- c(30, 60)
  plot_hlines <- NULL
  
  # Plot
  p <- ggplot(dt_edit, aes(x = value, y = key)) +
    # Content
    geom_bar(stat = "identity") + 
    geom_text(aes(label = value), 
              hjust = -.2,
              size = 2) + 

    # Faceting
    facet_grid(cols = vars(cond),
               scales = "free") +
    # Cosmetic
    labs(title = element_blank(),
         x     = element_blank(), 
         y     = element_blank()) + 
    theme(axis.text.x  = element_text(size = 5),
          axis.text.y  = element_text(size = 5),
          plot.margin  = unit(c(.05, .0, .0, .0), "cm"),
          # Background
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          # Facet Related
          strip.text = element_text(size = 8,
                                    face = "plain",
                                    margin = unit(c(.10, .10, .10, .10), "cm")) 
    ) + 
    # X Axis
    scale_x_continuous(breaks = plot_breaks,
                       labels = plot_labels) +
    geom_vline(xintercept = plot_vlines,
               size = .375,
               linetype = "dashed", color = "black") +
    coord_cartesian(xlim = plot_xlim) + 
    # Horizontal lines
    geom_hline(yintercept = plot_hlines, 
               size = .25,
               color = "gray")
  p
  return(p)
}

# ggplot facet_warp modify ------------------------------------------------
# source: https://fishandwhistle.net/post/2018/modifying-facet-scales-in-ggplot2/

# Inner function
scale_override <- function(which, scale) {
  if(!is.numeric(which) || (length(which) != 1) || (which %% 1 != 0)) {
    stop("which must be an integer of length 1")
  }
  
  if(is.null(scale$aesthetics) || !any(c("x", "y") %in% scale$aesthetics)) {
    stop("scale must be an x or y position scale")
  }
  
  structure(list(which = which, scale = scale), class = "scale_override")
}

## FacetWrap Version ##

CustomFacetWrap <- ggproto(
  "CustomFacetWrap", FacetWrap,
  init_scales = function(self, layout, x_scale = NULL, y_scale = NULL, params) {
    # make the initial x, y scales list
    scales <- ggproto_parent(FacetWrap, self)$init_scales(layout, x_scale, y_scale, params)
    
    if(is.null(params$scale_overrides)) return(scales)
    
    max_scale_x <- length(scales$x)
    max_scale_y <- length(scales$y)
    
    # ... do some modification of the scales$x and scales$y here based on params$scale_overrides
    for(scale_override in params$scale_overrides) {
      which <- scale_override$which
      scale <- scale_override$scale
      
      if("x" %in% scale$aesthetics) {
        if(!is.null(scales$x)) {
          if(which < 0 || which > max_scale_x) stop("Invalid index of x scale: ", which)
          scales$x[[which]] <- scale$clone()
        }
      } else if("y" %in% scale$aesthetics) {
        if(!is.null(scales$y)) {
          if(which < 0 || which > max_scale_y) stop("Invalid index of y scale: ", which)
          scales$y[[which]] <- scale$clone()
        }
      } else {
        stop("Invalid scale")
      }
    }
    
    # return scales
    scales
  }
)

facet_wrap_custom <- function(..., scale_overrides = NULL) {
  # take advantage of the sanitizing that happens in facet_wrap
  facet_super <- facet_wrap(...)
  
  # sanitize scale overrides
  if(inherits(scale_overrides, "scale_override")) {
    scale_overrides <- list(scale_overrides)
  } else if(!is.list(scale_overrides) || 
            !all(vapply(scale_overrides, inherits, "scale_override", FUN.VALUE = logical(1)))) {
    stop("scale_overrides must be a scale_override object or a list of scale_override objects")
  }
  
  facet_super$params$scale_overrides <- scale_overrides
  
  ggproto(NULL, CustomFacetWrap,
          shrink = facet_super$shrink,
          params = facet_super$params
  )
}

# EXAMPLE USE
# set.seed(123)
# test_large_values_data <- 1:9 %>%
#   set_names() %>%
#   map_dfr(~data.frame(
#     x_variable = runif(1, 100, 20000) + runif(10, -100, 100),
#     y_variable = runif(10, 0, 1)
#   ), .id = "facet_name")
# 
# p_annoying_x_scale <- ggplot(test_large_values_data, aes(x_variable, y_variable)) +
#   geom_point() +
#   facet_wrap(~facet_name, scales = "free", ncol = 4)
# 
# p_annoying_x_scale +
#   facet_wrap_custom(~facet_name, scales = "free", ncol = 4, scale_overrides = list(
#     scale_override(1, scale_x_continuous(breaks = c(5750, 5900))),
#     scale_override(6, scale_x_continuous(breaks = c(17800, 17900)))
#   ))

## Grid Version ##

CustomFacetGrid <- ggproto(
  "CustomFacetGrid", FacetGrid,
  init_scales = function(self, layout, x_scale = NULL, y_scale = NULL, params) {
    # make the initial x, y scales list
    scales <- ggproto_parent(FacetGrid, self)$init_scales(layout, x_scale, y_scale, params)
    
    if(is.null(params$scale_overrides)) return(scales)
    
    max_scale_x <- length(scales$x)
    max_scale_y <- length(scales$y)
    
    # ... do some modification of the scales$x and scales$y here based on params$scale_overrides
    for(scale_override in params$scale_overrides) {
      which <- scale_override$which
      scale <- scale_override$scale
      
      if("x" %in% scale$aesthetics) {
        if(!is.null(scales$x)) {
          if(which < 0 || which > max_scale_x) stop("Invalid index of x scale: ", which)
          scales$x[[which]] <- scale$clone()
        }
      } else if("y" %in% scale$aesthetics) {
        if(!is.null(scales$y)) {
          if(which < 0 || which > max_scale_y) stop("Invalid index of y scale: ", which)
          scales$y[[which]] <- scale$clone()
        }
      } else {
        stop("Invalid scale")
      }
    }
    
    # return scales
    scales
  }
)

facet_grid_custom <- function(..., scale_overrides = NULL) {
  # take advantage of the sanitizing that happens in facet_wrap
  facet_super <- facet_grid(...)
  
  # sanitize scale overrides
  if(inherits(scale_overrides, "scale_override")) {
    scale_overrides <- list(scale_overrides)
  } else if(!is.list(scale_overrides) || 
            !all(vapply(scale_overrides, inherits, "scale_override", FUN.VALUE = logical(1)))) {
    stop("scale_overrides must be a scale_override object or a list of scale_override objects")
  }
  
  facet_super$params$scale_overrides <- scale_overrides
  
  ggproto(NULL, CustomFacetGrid,
          shrink = facet_super$shrink,
          params = facet_super$params
  )
}

# EXAMPLE USE
# p_annoying_x_scale <- ggplot(test_large_values_data, aes(x_variable, y_variable)) +
#   geom_point() +
#   facet_grid(cols = vars(facet_name), scales = "free")
#
# p_annoying_x_scale <- ggplot(test_large_values_data, aes(x_variable, y_variable)) +
#   geom_point() +
#   facet_grid_custom(cols = vars(facet_name), scales = "free", scale_overrides = list(
#     scale_override(1, scale_x_continuous(breaks = c(5750, 5900))),
#     scale_override(6, scale_x_continuous(breaks = c(17800, 17900)))
#   ))
# Add a distinction
# head(test_large_values_data)
# dim(test_large_values_data)
# test_large_values_data$group <- 1:3
#
# p_annoying_x_scale <- ggplot(test_large_values_data, aes(x_variable, y_variable)) +
#   geom_point() +
#   facet_grid(rows = vars(group),
#              cols = vars(facet_name),
#              scales = "free")
# p_annoying_x_scale <- ggplot(test_large_values_data, aes(x_variable, y_variable)) +
#   geom_point() +
#   facet_grid_custom(rows = vars(group),
#                     cols = vars(facet_name),
#                     scales = "free", scale_overrides = list(
#     scale_override(1, scale_x_continuous(breaks = c(5750, 5900))),
#     scale_override(6, scale_x_continuous(breaks = c(17800, 17900)))
#   ))

plotwise <- function(res,
                     model,
                     parPlot,
                     meth_compare,
                     item_group){
  # Inputs
  # res = output
  # model = "sem"
  # parPlot = list(Means = 1:6,
  #                Variances = 7:12,
  #                Covariances = 13:27)
  # item_group = c(1:3) # items in a group recieving miss values
  # meth_compare = c("DURR_la","IURR_la", "blasso",
  #                  "MI_PCA", "MI_CART" ,"MI_RF",
  #                  "MI_OP",
  #                  "CC", "GS")

  ## Body
  # Counts
  n_conds <- nrow(res$conds)
  n_parms <- length(parPlot)

  # Extract Info from analysis model for both bias and CI coverage
  dt_bias = lapply(1:length(res[[model]]),
                   function(x) data.frame( res[[model]][[x]]$bias_per))[1:4]
  dt_cico = lapply(1:length(res[[model]]),
                   function(x) data.frame( res[[model]][[x]]$ci_cov))[1:4]

  # Isolate parameters and methods for comparison
  dt_bias <- lapply(parPlot, function(x){
    lapply(dt_bias, function(d){
      d[x, meth_compare]
    })
  })

  dt_cico <- lapply(parPlot, function(x){
    lapply(dt_cico, function(d){
      d[x, meth_compare]
    })
  })

  # Define Condition and Parameters names
  label_cond <- paste0("p = ",  res$conds$p, " ",
                       "pm = ", res$conds$pm)
  label_parm <- paste0(names(parPlot),
                       " (", sapply(parPlot, length), ")")

  # Create a data matrix structure for facet_grid
  store <- data.frame(parm = NA,
                      cond = NA,
                      methods = NA,
                      analysis = NA,
                      Min. = NA,
                      Mean = NA,
                      Max. = NA,
                      mean_i1 = NA,
                      mean_i2 = NA,
                      mean_i12 = NA
  )[0, ]

  # Populate matrix for facet_grid
  for (p in 1:n_parms) {
    # p <- 3
    for(cc in 1:n_conds){
      # cc <- 1
      # Extract Results to be processed
      results_bias <- abs(dt_bias[[p]][[cc]])
      results_cico <- abs(dt_cico[[p]][[cc]])

      # Create groups of parameters
      if(p != 3){
        par_group <- rep(1:2, each = length(parPlot[[p]])/2)
      }
      if(p == 3){
        # Define the unique covariances
        x <- t(combn(1:length(parPlot$Means), 2))
        par_group <- rep(NA, nrow(x))
        # Define which unique covariances go where
        par_group[which(x[, 1] %in% item_group & x[, 2] %in% item_group)] <- 1
        par_group[which(x[, 1] %in% item_group & !x[, 2] %in% item_group)] <- 2
        par_group[which(!x[, 1] %in% item_group & !x[, 2] %in% item_group)] <- 3
      }

      # Obtain Min, Mean, Max
      out_bias <- t(sapply(results_bias, summary)[c("Min.", "Mean", "Max."), ])
      out_cico <- t(sapply(results_cico, summary)[c("Min.", "Mean", "Max."), ])

      # Group Means
      out_bias <- cbind(out_bias, t(sapply(results_bias, tapply, par_group, mean)))
      out_cico <- cbind(out_cico, t(sapply(results_cico, tapply, par_group, mean)))

      # Stack results
      results <- rbind(out_bias, out_cico)

      # Correct names and dimensions
      rownames(results) <- sub("_", "-", rownames(results))
      if(ncol(results) < 6){
        results <- cbind(results, NA)
      }

      # Attach Methods as factor w/ precise levels order
      out <- data.frame(
        parm = label_parm[p],
        cond = label_cond[cc],
        methods = factor(rownames(results),
                         levels = rownames(results),
                         labels = rownames(results)),
        analysis = rep(c("Bias", "CI coverage"), each = nrow(out_bias)),
        results
      )
      colnames(out) <- colnames(store)

      # Attach Condition Name
      store <- rbind(store, out)
    }
  }

  # Reshape
  store_melt <- reshape2::melt(store,
                               id.var = c("parm", "cond", "methods", "analysis"))
  levels(store_melt$variable) <- c("Min", "Mean", "Max",
                                   "Item set 1", "Item set 2",
                                   "Item set 1 and 2")
  # Return data.frame Output
  return(store_melt)
}