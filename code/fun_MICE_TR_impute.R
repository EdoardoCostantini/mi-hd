### Title:    Imputation according to optimal oracle imputation model
### Project:  Imputing High Dimensional Data
### Author:   Anonymized for peer review
### Created:  2020-05-19

impute_MICE_OP <- function(Z, O, cond, perform = TRUE, parms = parms){
  
  ## Input: 
  # @Z: dataset w/ missing values
  # @AS_size a number indicating the size of the Active Set
  # @chains: number of imputation chains
  # @iters: itnerations and number
  # @parms: the initialization object parms
  
  ## Example inputs
  # Z = Xy_mis # or
  # Z = Xy_MIOP[, CIDX_MOP]
  # O <- as.data.frame(!is.na(Z))            # matrix index of observed values
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
    
    # Define predictor matrix for MI TRUE with best active set
    p_imp_id <- which(colSums(O) != nrow(O))
    predMat <- matrix(rep(0, ncol(Z)^2), ncol = ncol(Z), 
                      dimnames = list(colnames(Z), colnames(Z)))
    var_index <- colnames(predMat) %in% parms$S_all # flexibility interaction!
    predMat[p_imp_id, var_index] <- 1
    
    # Other variables needing imputation
    # Separate variables target prespecified variable selection for Imp Mods
    target <- which(colnames(Z) %in% parms$z_m_id)
    if(sum(is.na(Z[, -target])) > 0){
      pMat     <- quickpred(Z, mincor = .3)
      
      # Define Single Imputed data as the auxiliary set
      predMat[-target, ] <- pMat[-target, ]
    }

    # Define methods
    methods <- rep("norm", ncol(Z))
    vartype <- sapply(Z, class)
    methods[vartype != "numeric"] <- "pmm"
    
    # For Experiment 4: optimal predictors are found with mice function
    if(parms$exp == 4){
      predMat <- quickpred(Z,
                           mincor = .3,
                           include = parms$S_all)
      methods <- "pmm"
    }
    
    # Impute
    imp_MITR_mids <- mice::mice(Z, 
                                predictorMatrix = predMat,
                                m = parms$mice_ndt,
                                maxit = parms$mice_iters,
                                printFlag       = FALSE,
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
