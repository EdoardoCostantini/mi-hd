### Title:    Imputing High Dimensional Data: cart imputation
### Author:   Edoardo Costantini
### Created:  2020-05-19

impute_MICE_CART <- function(Z, parms){
  
  ## Input: 
  # Z = Xy_mis # dataset w/ missing values, 
  
  ## output: 
  # - a list of chains imputed datasets at iteration iters, 
  # - per variable list of imputed values
  # - imputation run time
  
  ## body:
  if(parms$meth_sel$MI_CART == TRUE){
    start.time <- Sys.time()
    
    imp_MI_CART_mids <- mice::mice(Z, 
                                   m = parms$mice_ndt, 
                                   maxit = parms$mice_iters,
                                   meth = 'cart', 
                                   minbucket = 5)
    
    end.time <- Sys.time()
    
    imp_MICART_dats <- mice::complete(imp_MI_CART_mids, "all")
    imp_MICART_imps <- imp_MI_CART_mids$imp
    imp_MICART_time <- difftime(end.time, start.time, units = "mins")
    
    return(list(dats = imp_MICART_dats,
                imps = imp_MICART_imps,
                time = imp_MICART_time))
  } else {
    return(list(dats = NULL,
                imps = NULL,
                time = NULL))
  }
}

impute_MICE_CARTbb <- function(Z, parms){
  
  ## Input: 
  # @Z: dataset w/ missing values, 
  
  ## output: 
  # - a list of chains imputed datasets at iteration iters, 
  # - per variable list of imputed values
  # - imputation run time
  
  ## body:
  if(parms$meth_sel$MI_CTBB == TRUE){
    start.time <- Sys.time()
    
    imp_MI_CART_mids <- miceImpHDv::mice(Z, 
                                         m = parms$ndt, 
                                         maxit = parms$iters,
                                         meth = 'cart.bb', 
                                         minbucket = 5)
    
    end.time <- Sys.time()
    
    imp_MICART_dats <- mice::complete(imp_MI_CART_mids, "all")
    imp_MICART_imps <- imp_MI_CART_mids$imp
    imp_MICART_time <- difftime(end.time, start.time, units = "mins")
    
    return(list(dats = imp_MICART_dats,
                imps = imp_MICART_imps,
                time = imp_MICART_time))
  } else {
    return(list(dats = NULL,
                imps = NULL,
                time = NULL))
  }
}
