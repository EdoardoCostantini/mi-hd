### Title:    Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2021-02-27
### Description: This is the PcAux version for the project of the 
###              MI-PCA.

impute_PcAux <- function(Z, O, cond, DA = FALSE, parms = parms){
  
  ## Input: 
  # @Z: dataset w/ missing values, 
  # @AS_size a number indicating the size of the Active Set
  # @chains: number of imputation chains, 
  # @iters: itnerations and number
  # @parms: the initialization object parms
  
  ## output: 
  # - a list of chains imputed datasets at iteration iters
  # - per variable list of imputed values
  # - imputation run time
  
  # For internals
  # DA = cond$int_da
  # Z = Xy_mis
  # O = as.data.frame(!is.na(Z)) # matrix index of observed values
  # DA = FALSE
  
  
  ## body:
  if(parms$meth_sel$MI_PCA == TRUE){

    tryCatch({
      
      start.time <- Sys.time()

      # Separate variables target of imputation from auxiliary variables 
      target <- which(colnames(Z) %in% parms$z_m_id)

      # Prepare Data Object
      cleanAux <- prepData(rawData = Z,
                           simMode = TRUE,
                           verbose = 1)
      
      # Extracting PCs
      print("PCA Impute: Extracting PCs")
      pcAuxDat <- createPcAux(pcAuxData    = cleanAux,
                              interactType = 0,    # no interaction terms
                              maxPolyPow   = 1L,   # no polynomials
                              simMode      = TRUE, # suppress automiatic data checks
                              pcaMemLevel  = 1L,
                              nComps       = c(0.5, 0),
                              control = list(
                                miceIters = 10
                              )
      )
      
      # Store Components
      Z_pcs <- pcAuxDat$pcAux$lin
      
      # Merge to original data
      Z_imp_input <- cbind(Z, Z_pcs)
      
      # Define Imputation Model Predictor Matrix
      predMat <- matrix(rep(0, (ncol(Z_imp_input)^2)), 
                        ncol = ncol(Z_imp_input), 
                        dimnames = list(colnames(Z_imp_input), 
                                        colnames(Z_imp_input)
                        )
      )
      
      PC_colIndex <- (ncol(Z)+1):(ncol(Z)+ncol(Z_pcs))
      predMat[, PC_colIndex] <- 1
      
      # Perform Imputation
      print("PCA Impute: Performing MI")
      imp_PCA_mids <- mice::mice(Z_imp_input,
                                 m      = parms$mice_ndt,
                                 maxit  = parms$mice_iters,
                                 predictorMatrix = predMat,
                                 printFlag = TRUE,
                                 method = "pmm")
      imp_PCA_mids$predictorMatrix[target, c(1:10, 95:105, 500:527)]
      
      # Store results
      print("PCA Impute: Storing Results")
      imp_PCA_PC_dats <- mice::complete(imp_PCA_mids, "all")
      
      # Fill into orignal datasets the imputations (get rid of PCs)
      imp_PCA_dats <- lapply(imp_PCA_PC_dats, function(x){
        x <- x[, -PC_colIndex]
        return(x)
      })

      imp_PCA_imps <- imp_PCA_mids$imp
      imp_PCA_time <- difftime(end.time, start.time, units = "mins")
      
      return(list(dats = imp_PCA_dats,
                  imps = imp_PCA_imps,
                  time = imp_PCA_time,
                  mids = imp_PCA_mids,
                  dtIN = Z_imp_input))
      
      ### END TRYCATCH EXPRESSION
    }, error = function(e){
      err <- paste0("Original Error: ", e)
      print(err)
      return(list(dats = NULL,
                  imps = NULL,
                  time = NULL,
                  mids = NULL,
                  dtIN = NULL)
      )
    }
    )
    
  } else {
    return(list(dats = NULL,
                imps = NULL,
                time = NULL,
                mids = NULL,
                dtIN = NULL))
  }
}
