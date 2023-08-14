### Title:    Imputation model defined based on quickpred()
### Project:  Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2022-01-28

impute_MICE_qp <- function(Z, perform = TRUE, parms = parms,
                           ridge = 0, eps = 0, threshold = 1){
  
  ## Input: 
  # @Z: dataset w/ missing values
  # @parms: the initialization object parms
  
  ## Example inputs
  # Z = Xy_mis # or
  # parms = parms
  
  ## output: 
  # - a list of chains imputed datasets at iteration iters
  # - per variable list of imputed values
  # - imputation run time
  
  ## body:
  if(perform == TRUE){
    
    tryCatch({
      
    start.time <- Sys.time()
    
    # Define predictor matrix based on quickpred
    pMat <- quickpred(Z)

    # Define methods
    methods <- rep("norm", ncol(Z))
    vartype <- sapply(Z, class)
    methods[vartype != "numeric"] <- "pmm"

    # Impute
    imp_MIqp_mids <- mice::mice(Z,
                                predictorMatrix = pMat,
                                m = parms$mice_ndt,
                                maxit = parms$mice_iters,
                                printFlag = TRUE,
                                ridge = ridge,
                                eps = eps,
                                threshold = threshold,
                                method = methods)
    end.time <- Sys.time()
    
    # Store results
    imp_MIqp_dats <- mice::complete(imp_MIqp_mids, "all")
    imp_MIqp_imps <- imp_MIqp_mids$imp
    imp_MIqp_time <- difftime(end.time, start.time, units = "mins")
    
    return(list(dats = imp_MIqp_dats,
                imps = imp_MIqp_imps,
                time = imp_MIqp_time,
                mids = imp_MIqp_mids))
    
    
    }, error = function(e){
      err <- paste0("Original Error: ", e)
      print(err)
      return(list(dats = NULL,
                  imps = NULL,
                  time = NULL,
                  mids = NULL)
      )
    }
    )

  } else {
    return(list(dats = NULL,
                imps = NULL,
                time = NULL,
                mids = NULL))
  }
}
