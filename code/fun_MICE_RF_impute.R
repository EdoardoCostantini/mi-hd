### Title:    Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-05-19

impute_MICE_RF <- function(Z, parms){
  
  ## Input: 
  # @Z: dataset w/ missing values, 
  
  ## output: 
  # - a list of chains imputed datasets at iteration iters, 
  # - per variable list of imputed values
  # - imputation run time
  
  ## body:
  if(parms$meth_sel$MI_RF == TRUE){
    
    tryCatch({
      
      start.time <- Sys.time()
      imp_MI_RF_mids <- mice::mice(Z, 
                                   m = parms$mice_ndt,
                                   maxit = parms$mice_iters,
                                   meth = "rf", 
                                   ntree = parms$rfntree)
      end.time <- Sys.time()
      
      imp_MIRF_dats <- mice::complete(imp_MI_RF_mids, "all")
      imp_MIRF_imps <- imp_MI_RF_mids$imp
      imp_MIRF_time <- difftime(end.time, start.time, units = "mins")
      
      return(list(dats = imp_MIRF_dats,
                  imps = imp_MIRF_imps,
                  time = imp_MIRF_time))
    ### END TRYCATCH EXPRESSION
    }, error = function(e){
      err <- paste0("Original Error: ", e)
      print(err)
      return(list(dats = NULL,
                  imps = NULL,
                  time = NULL)
      )
    }
    )
    
  } else {
    return(list(dats = NULL,
                imps = NULL,
                time = NULL))
  }
}
