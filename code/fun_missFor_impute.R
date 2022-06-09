### Title:    Imputing High Dimensional Data
### Author:   Anonymized for peer review
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
    
    suppressWarnings(
      # In exp4, where ordered items w/ 4 categories need to be imputed
      # having the item as numeric passed to the randomForest function
      # produces a warning asking whether you are sure you want to use
      # regression. To use a decision tree you would need to pass the
      # y_obs argument as a factor. However, for consistency we stick
      # to the numeric treatment and suppress the warnings that come
      # out of it
      imp_missForest <- missForest(Z,
                                   maxiter = parms$missFor_maxiter,
                                   ntree = parms$missFor_ntree,
                                   # maxiter = 20,
                                   # ntree = 100,
                                   parallelize = "no")
    )
    
    end.time <- Sys.time()
    
    imp_missFor_dats <- imp_missForest$ximp
    
    return(list(dats = imp_missFor_dats))
    
  } else {
    return(list(dats = NULL))
  }
}
