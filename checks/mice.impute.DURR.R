### Title:    Adding my own imputation function
### Author:   Anonymized for peer review
### Created:  2020-06-08
### Notes:    Create a version of mice.impute.method that corresponds to your 
###           DURR proposition in DengEtAl 2016

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

mice.impute.DURR <- function (y, ry, x, wy = NULL, ...)
{
  if (is.null(wy))
    wy <- !ry
  
  Zm <- cbind(y, x)
  
  idx_bs   <- sample(1:nrow(Zm), nrow(Zm), replace = TRUE)
  Zm_bs <- Zm[idx_bs, ]
  
  # Select data
  y_obs_bs <- Zm_bs[ry[idx_bs], 1]  # z_j_obs
  y_mis_bs <- Zm_bs[!ry[idx_bs], 1] # zm_mj [useless]
  X_obs_bs <- Zm_bs[ry[idx_bs], -1] # Wm_j_obs
  X_mis <- Zm[!ry, -1] # Wm_mj
  
  glmfam <- detect_family(y)
  reg_type <- "lasso"
  
  # 2. Fit regularized regression on bootstraped observed data
  ## Lasso
  if(reg_type == "lasso"){
    regu.mod <- rr_est_lasso(X = X_obs_bs, y = y_obs_bs, 
                             parms = parms, fam = glmfam)
  }
  
  ## Elastic net
  if(reg_type == "el"){
    regu.mod <- rr_est_elanet(X = X_obs_bs, y = y_obs_bs, 
                              parms = parms, fam = glmfam)
  }
  
  # 3. Predict zm_j based on onrignal data (not bootsrap)
  # (i.e. obtain imputations (imps))
  if(glmfam == "gaussian"){
    zm_j <- imp_gaus_DURR(regu.mod, X_obs_bs, y_obs_bs, X_mis, parms)
  }
  if(glmfam == "binomial"){
    zm_j <- imp_dich_DURR(regu.mod, X_obs_bs, y_obs_bs, X_mis, parms)
  }
  if(glmfam == "multinomial"){
    zm_j <- imp_multi_DURR(regu.mod, X_obs_bs, y_obs_bs, X_mis, parms)
  }
  
  return(zm_j)
}

