### Title:    Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-05-19

impute_MICE_TR <- function(Z, cond, perform = TRUE, parms = parms){
  
  ## Input: 
  # @Z: dataset w/ missing values, 
  # @AS_size a number indicating the size of the Active Set
  # @chains: number of imputation chains, 
  # @iters: itnerations and number
  # @parms: the initialization object parms
  
  ## Example inputs
  # Z = Xy_mis
  # cond = conds[1,]
  # parms = parms
  
  ## output: 
  # - a list of chains imputed datasets at iteration iters
  # - per variable list of imputed values
  # - imputation run time
  
  ## body:
  if(perform == TRUE){
    # Select true active set 
    S <- parms$S_all
    MI_ture_pred <- c(S)
    
    # Define predictor matrix for MI TRUE
    predMat <- matrix(rep(0, ncol(Z)^2), ncol = ncol(Z), 
                      dimnames = list(colnames(Z), colnames(Z)))
    predMat[1:length(parms$z_m_id), MI_ture_pred] <- 1
    
    # Define methods
    methods <- rep("norm", ncol(Z))
    vartype <- sapply(Z, class)
    methods[vartype != "numeric"] <- "pmm"
    
    # Impute
    start.time <- Sys.time()
    
    imp_MITR_mids <- mice::mice(Z, 
                                predictorMatrix = predMat,
                                m = parms$mice_ndt,
                                maxit = parms$mice_iters,
                                ridge = 1e-5,
                                method = methods)
    
    end.time <- Sys.time()
    
    # Store results
    imp_MITR_dats <- mice::complete(imp_MITR_mids, "all")
    imp_MITR_imps <- imp_MITR_mids$imp
    imp_MITR_time <- difftime(end.time, start.time, units = "mins")
    
    return(list(dats = imp_MITR_dats,
                imps = imp_MITR_imps,
                time = imp_MITR_time,
                mids = imp_MITR_mids))
  } else {
    return(list(dats = NULL,
                imps = NULL,
                time = NULL,
                mids = NULL))
  }
}
