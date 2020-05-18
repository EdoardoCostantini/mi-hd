### Title:    imputeHD-comp impute w/ Regu Freq Regr DURR
### Author:   Edoardo Costantini
### Created:  2020-05-08

impute_MICE_TR <- function(Z, AS_size, parms = parms){
  
  ## Input: 
  # @Z: dataset w/ missing values, 
  # @AS_size a number indicating the size of the Active Set
  # @chains: number of imputation chains, 
  # @iters: itnerations and number
  # @parms: the initialization object parms
  
  ## output: 
  # - a list of chains imputed datasets at iteration iters
  # - per variable list of imputed values
  # - imputation run time
  
  ## body:
  # Select true active set 
  S <- parms$S_all[[ which(paste0("q", AS_size) == names(parms$S_all)) ]]
  MI_ture_pred <- c((S), ncol(Z))
  
  # Define predictor matrix for MI TRUE
  predMat <- matrix(rep(0, ncol(Z)^2), ncol = ncol(Z), 
                    dimnames = list(colnames(Z), colnames(Z)))
  predMat[parms$z_m_id, MI_ture_pred] <- 1
  # predMat[c("z1", "z2", "z3"), MI_ture_pred] <- 1
  
  # Impute
  start.time <- Sys.time()
  
  imp_MITR_mids <- mice::mice(Z, 
                              predictorMatrix = predMat,
                              m = parms$ndt,
                              maxit = parms$iters,
                              ridge = 1e-5,
                              method = "norm")
  
  end.time <- Sys.time()
  
  # Store results
  imp_MITR_dats <- mice::complete(imp_MITR_mids, "all")
  imp_MITR_imps <- imp_MITR_mids$imp
  imp_MITR_time <- difftime(end.time, start.time, units = "mins")
  
  return(list(dats = imp_MITR_dats,
              imps = imp_MITR_imps,
              time = imp_MITR_time,
              mids = imp_MITR_mids))
}
