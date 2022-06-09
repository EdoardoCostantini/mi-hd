### Title:    Elastic Net imputation time and cross-validation
### Project:  Imputing High Dimensional Data
### Author:   Anonymized for peer review
### Created:  2020-05-19

rm(list=ls())
source("./init.R")
source("./functions.R")

# Missing data pm ---------------------------------------------------------
# Does the proportion of missing data imposed work?

Xy <- genDt_mvn(conds[1, ], parms)
Xy_mis <- imposeMiss(Xy, parms, conds[1, ])
Xy_mis <- cbind(Xy_mis[parms$z_m_id], Xy_mis[-parms$z_m_id])

O <- !is.na(Xy_mis)

# -------------------------------------------------------------------------
# Check how much time does one iteration take
parms$iters <- 5
parms$burnin_imp <- 0
parms$ndt <- 5

DURR_time <- system.time(
  imp_DURR_el <- impute_DURR(Z = Xy_mis,
                             O = as.data.frame(O),
                             reg_type = "el",
                             perform = TRUE,
                             parms = parms)
)

DURR_time/parms$iters


# -------------------------------------------------------------------------
# Check different types of corssvalidation.

## GRID ##


## RANDOM ##
  
  Z = Xy_mis
  O <- as.data.frame(!is.na(Xy_mis))            # matrix index of observed values
  reg_type = c("el", "lasso")[1] # imputation model penality type
  
  p  <- ncol(Z) # number of variables [indexed with j]
  
  p_imp <- length(parms$z_m_id) # variables with missing values
  r     <- colSums(O[1:p_imp]) # variable-wise count of observed values
  nr    <- colSums(!O[1:p_imp]) # variable-wise count of observed values
  
  start.time <- Sys.time()
  
  Zm <- init_dt_i(Z, missing_type(Z)) # reinitialize data
  
  j <- 1
  # 1. Generate bootstrap data set (for each j!)
  # Generate bootstrap sample
  idx_bs   <- sample(1:nrow(Zm), nrow(Zm), replace = TRUE)
  Zm_bs <- Zm[idx_bs, ]
  
  # Select data
  y_obs_bs <- Zm_bs[O[idx_bs, j], j]  # z_j_obs
  y_mis_bs <- Zm_bs[!O[idx_bs, j], j] # zm_mj [useless]
  X_obs_bs <- as.matrix(Zm_bs[O[idx_bs, j], -j]) # Wm_j_obs
  X_obs_bs <- apply(X_obs_bs, 2, as.numeric) # makes dicho numbers
  X_mis <- Wm_mj    <- as.matrix(Zm[!O[, j], -j]) # Wm_mj
  X_mis <- apply(X_mis, 2, as.numeric)
  
  glmfam <- detect_family(Zm[, j])
  
  # 2. Fit regularized regression on bootstraped observed data
  ## Elastic net
  if(reg_type == "el"){
    regu.mod <- rr_est_elanet(X = X_obs_bs, y = y_obs_bs, 
                              parms = parms, fam = glmfam)
  }
  
  # Set training control
  train_control <- trainControl(method = "cv",
                                number = 10, # 10-fold cross validation 
                                selectionFunction = "best",
                                verboseIter = FALSE)
  
  # CV
  el_cv <- train(y_obs_bs ~ .,
                 data = data.frame(y_obs_bs, X_obs_bs), 
                 method = "glmnet",
                 family = glmfam, 
                 type.multinomial = "grouped",
                 metric = "RMSE",
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
  
  # 3. Predict zm_j based on onrignal data (not bootsrap)
  # (i.e. obtain imputations (imps))
  
  zm_j <- imp_gaus_DURR(regu.mod, X_obs_bs, y_obs_bs, X_mis, parms)
  
  end.time <- Sys.time()
