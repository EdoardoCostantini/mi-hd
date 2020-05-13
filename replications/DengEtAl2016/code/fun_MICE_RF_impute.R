### Title:    imputeHD-comp impute w/ Regu Freq Regr DURR
### Author:   Edoardo Costantini
### Created:  2020-05-08

impute_MICE_RF <- function(Xy_mis, chains=5, iters=5, rfntree = 10){
  
  ## Input: 
  # @Xy_mis: dataset w/ missing values, 
  # @chains: number of imputation chains, 
  # @iters: itnerations and number
  # @rfntree: number of trees for random forest
  
  ## output: 
  # - a list of chains imputed datasets at iteration iters, 
  # - per variable list of imputed values
  # - imputation run time
  
  ## body:
  start.time <- Sys.time()
  imp_MI_RF_mids <- mice::mice(Xy_mis, 
                               m = chains,
                               maxit = iters,
                               meth = "rf", 
                               ntree = rfntree)
  end.time <- Sys.time()
  
  imp_MIRF_dats <- mice::complete(imp_MI_RF_mids, "all")
  imp_MIRF_imps <- imp_MI_RF_mids$imp
  imp_MIRF_time <- difftime(end.time, start.time, units = "mins")
  
  return(list(dats = imp_MIRF_dats,
              imps = imp_MIRF_imps,
              time = imp_MIRF_time))
}
