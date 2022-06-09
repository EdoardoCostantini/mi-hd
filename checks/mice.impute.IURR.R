### Title:    mice.impute.IURR
### Author:   Anonymized for peer review
### Created:  2020-06-09
### Notes:    Create a version of mice.impute.method that corresponds to your 
###           IURR proposition in DengEtAl 2016

# # Copy an existing imputation function
# # mice.impute.norm
# function (y, ry, x, wy = NULL, ...)
# {
#   if (is.null(wy))
#     wy <- !ry
#   x <- cbind(1, as.matrix(x))
#   parm <- .norm.draw(y, ry, x)
#   return(x[wy, ] %*% parm$beta + rnorm(sum(wy)) * parm$sigma)
# }

mice.impute.IURR <- function (y, ry, x, wy = NULL, ...)
{
  if (is.null(wy))
    wy <- !ry
  
  Zm <- cbind(y, x)
  
  # Select data
  y_obs <- Zm[ry, 1]  # z_j_obs
  y_mis <- Zm[!ry, 1] # zm_mj [useless]
  X_obs <- Zm[ry, -1] # Wm_j_obs
  X_mis <- Zm[!ry, -1] # Wm_mj
  
  glmfam <- detect_family(y)
  reg_type <- "lasso"
  
  # 2. Fit regularized regression on bootstraped observed data
  IURR_fwd_cond <- FALSE
  while (IURR_fwd_cond == FALSE) {
    
    # Lasso
    if(reg_type == "lasso"){
      regu.mod <- rr_est_lasso(X = X_obs, y = y_obs, 
                               parms = parms, fam = glmfam)
    }
    
    # Elastic net
    if(reg_type == "el"){
      regu.mod <- rr_est_elanet(X = X_obs, y = y_obs, 
                                parms = parms, fam = glmfam)
    }
    
    IURR_fwd_cond <- length(y_obs) > sum(coef(regu.mod) != 0)
    
  }
  
  # 3. Predict zm_j (i.e. obtain imputations (imps))
  if(glmfam == "gaussian"){
    zm_j <- imp_gaus_IURR(model = regu.mod,
                          X_tr = X_obs, y_tr = y_obs,
                          X_te = X_mis, y_te = y_mis, 
                          parms = parms)
  }
  if(glmfam == "binomial"){
    zm_j <- imp_dich_DURR(model = regu.mod,
                          X_tr = X_obs, y_tr = y_obs,
                          X_te = X_mis, 
                          parms = parms)
  }
  if(glmfam == "multinomial"){
    zm_j <- imp_multi_DURR(regu.mod, X_obs_bs, y_obs_bs, X_mis, parms)
  }
  
  return(zm_j)
}

