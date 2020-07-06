### Title:    Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-05-19

impute_missFor <- function(Z, parms){
  
  ## Input: 
  # @Z: dataset w/ missing values, 
  
  ## output: 
  # - a list of chains imputed datasets at iteration iters, 
  # - per variable list of imputed values
  # - imputation run time
  
  ## body:
  if(parms$meth_sel$missFor == TRUE){
    start.time <- Sys.time()
    
    imp_missForest <- missForest(Z,
                                 maxiter = 20,
                                 ntree = 100,
                                 parallelize = "no")
    
    end.time <- Sys.time()
    
    imp_missFor_dats <- imp_missForest$ximp
    
    return(list(dats = imp_missFor_dats))
    
  } else {
    return(list(dats = NULL))
  }
}
