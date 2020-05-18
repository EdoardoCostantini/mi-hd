### Title:    imputeHD-comp impute w/ Regu Freq Regr DURR
### Author:   Edoardo Costantini
### Created:  2020-05-08

impute_MICE_RF <- function(Z, parms){
  
  ## Input: 
  # @Z: dataset w/ missing values, 
  
  ## output: 
  # - a list of chains imputed datasets at iteration iters, 
  # - per variable list of imputed values
  # - imputation run time
  
  ## body:
  start.time <- Sys.time()
  imp_MI_RF_mids <- mice::mice(Z, 
                               m = parms$ndt,
                               maxit = parms$iters,
                               meth = "rf", 
                               ntree = parms$rfntree)
  end.time <- Sys.time()
  
  imp_MIRF_dats <- mice::complete(imp_MI_RF_mids, "all")
  imp_MIRF_imps <- imp_MI_RF_mids$imp
  imp_MIRF_time <- difftime(end.time, start.time, units = "mins")
  
  return(list(dats = imp_MIRF_dats,
              imps = imp_MIRF_imps,
              time = imp_MIRF_time))
}
