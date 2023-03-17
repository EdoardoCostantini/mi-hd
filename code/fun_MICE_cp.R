### Title:    Imputation model using custom predictors
### Project:  Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2022-01-28

impute_MICE_cp <- function(Z, preds = c(), perform = TRUE, parms = parms,
                           ridge = 0, eps = 0, threshold = 1){
  
  ## Input: 
  # @Z: dataset w/ missing values
  # @parms: the initialization object parms
  
  ## Example inputs
  # Z = Xy_mis # or
  # preds = parms$z_m_id
  # cond = conds[2,]
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
    pMat <- make.predictorMatrix(Z)
    pMat[, -which(colnames(Z) %in% preds)] <- 0

    # Define methods
    methods <- rep("norm", ncol(Z))
    vartype <- sapply(Z, class)
    methods[vartype != "numeric"] <- "pmm"

    # Impute
    imp_MIcp_mids <- mice::mice(Z,
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
    imp_MIcp_dats <- mice::complete(imp_MIcp_mids, "all")
    imp_MIcp_imps <- imp_MIcp_mids$imp
    imp_MIcp_time <- difftime(end.time, start.time, units = "mins")
    
    return(list(dats = imp_MIcp_dats,
                imps = imp_MIcp_imps,
                time = imp_MIcp_time,
                mids = imp_MIcp_mids))
    
    
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
