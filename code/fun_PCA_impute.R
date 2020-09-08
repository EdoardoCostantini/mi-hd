### Title:    Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-05-19

impute_PCA <- function(Z, O, cond, DA = FALSE, parms = parms){
  
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
  # Z = Xy_input[, CIDX_all]
  # Z = Xy_mis
  # O = as.data.frame(!is.na(Z)) # matrix index of observed values
  # DA = cond$int_da
  
  ## body:
  if(parms$meth_sel$MI_PCA == TRUE){
    start.time <- Sys.time()
    
    # Data
    p <- ncol(Z) # number of variables [INDEX with j]
    Z_aux <- Z
    
    # Separate variables target of imputation from auxiliary variables 
    target <- which(colnames(Z) %in% parms$z_m_id)
    
    # Generate DA version of Z_aux if required
    if(DA == TRUE){
      # Create an augmented Z dataset (inlcuding ALL interaction)
      print("PCA Impute: Creating Two Way Interactions")
      twoWays <- computeInteract(scale(Z,
                                       center = TRUE,
                                       scale  = FALSE),
                                 idVars = colnames(Z),
                                 ordVars = NULL,
                                 nomVars = NULL,
                                 moderators = colnames(Z))
      Z_aux <- cbind(Z, twoWays)
    }
    
    if(sum(is.na(Z_aux[, -target])) > 0){
      # Fill in cases when ZDA has missing values
      print("PCA Impute: Filling in auxiliary NAs")
      
      pMat     <- quickpred(Z_aux, mincor = .3)
      ZDA_mids <- mice(Z_aux,
                       m               = 1, 
                       maxit           = parms$SI_iter,
                       predictorMatrix = pMat,
                       printFlag       = FALSE,
                       ridge           = cond$ridge,
                       method          = "pmm")
      Z_out <- complete(ZDA_mids)
      
      # Define Single Imputed data as the auxiliary set
      Z_aux[, -target] <- Z_out[, -target]
      
      # Define dataset for output (will be used by other imputation methods)
      Z_out <- Z_aux
      
    } else {
      Z_out <- Z_aux # Need it for output consistency
      
    }
    
    # Extract PCs
    pcaOut <- prcomp(Z_aux[, -target], scale = TRUE, retx = TRUE)
    
    ## Compute and store the cumulative proportion of variance explained by
    ## the component scores:
    rSquared <- cumsum(pcaOut$sdev^2) / sum(pcaOut$sdev^2)
    
    ## Extract the principal component scores:
    Z_pca <- pcaOut$x[, rSquared < parms$PCA_pcthresh]
    
    ## Define Imputation methods
    Z_input <- cbind(Z[, target], Z_pca)
      methods <- rep("norm", ncol(Z_input))
      vartype <- sapply(Z_input, class)
      methods[vartype != "numeric"] <- "pmm"
      
    ## Impute
    print("PCA Impute: Performing Multiple Imputation")
    imp_PCA_mids <- mice::mice(Z_input,
                               m      = parms$mice_ndt,
                               maxit  = parms$mice_iters,
                               printFlag = FALSE,
                               ridge  = cond$ridge,
                               method = methods)
    
    end.time <- Sys.time()
    
    # Store results
    print("PCA Impute: Storing Results")
    imp_PCA_PC_dats <- mice::complete(imp_PCA_mids, "all")
    
    # Fill into orignal datasets the imputations (get rid of PCs)
    imp_PCA_dats <- lapply(imp_PCA_PC_dats, function(x){
      df_temp <- Z
      df_temp[, parms$z_m_id] <- x[, parms$z_m_id]
      return(df_temp)
    })
    
    imp_PCA_imps <- imp_PCA_mids$imp[1:parms$zm_n]
    imp_PCA_time <- difftime(end.time, start.time, units = "mins")
    
    return(list(dats = imp_PCA_dats,
                imps = imp_PCA_imps,
                time = imp_PCA_time,
                mids = imp_PCA_mids,
                dtIN = Z_out))
  } else {
    return(list(dats = NULL,
                imps = NULL,
                time = NULL,
                mids = NULL,
                dtIN = NULL))
  }
}
