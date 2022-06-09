### Title:    Imputing High Dimensional Data
### Author:   Anonymized for peer review
### Created:  2020-05-19

impute_PCA <- function(Z, O, cond, parms = parms){
  
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
  # O = as.data.frame(!is.na(Z)) # matrix index of observed values
  # Z = Xy_mis
  # Z = Xy_SI
  # O = as.data.frame(!is.na(Z)) # matrix index of observed values
  
  ## body:
  if(parms$meth_sel$MI_PCA == TRUE){

    tryCatch({
      
      start.time <- Sys.time()
      
      # Data
      p <- ncol(Z) # number of variables [INDEX with j]
      Z_aux <- Z
      
      # Separate variables target of imputation from auxiliary variables 
      target <- which(colnames(Z) %in% parms$z_m_id)
      
      if(sum(is.na(Z_aux[, -target])) > 0){
        # Fill in cases when ZDA has missing values
        print("PCA Impute: Filling in auxiliary NAs")
        
        pMat     <- quickpred(Z_aux, mincor = .3)
        ZDA_mids <- mice(Z_aux,
                         m               = 1, 
                         maxit           = 100,
                         predictorMatrix = pMat,
                         # printFlag       = FALSE,
                         printFlag       = TRUE,
                         ridge           = cond$ridge,
                         method          = "pmm")
        
        Z_out <- complete(ZDA_mids)
        
        # Define Single Imputed data as the auxiliary set
        Z_aux[, -target] <- Z_out[, -target]
        
        # Define dataset for output
        Z_out <- Z_aux

      } else {
        Z_out <- Z_aux # Need it for output consistency
        
      }
      # If there are categorical variables, then use FAMD
      if("factor" %in% sapply(Z_aux, class)){
        
        res.famd <- FAMD(Z_aux[, -target], graph = FALSE, ncp = ncol(Z_aux))
        Z_pca <- res.famd$ind$coord[, res.famd$eig[, 3] <= parms$PCA_pcthresh*100]
        
      } else {
      # If there are continuous and dichotomous variables, then use prcomp
        # Create Model matrix for PCA extraction
        Z_PC <- model.matrix(~ ., Z_aux[, -target])[, -1]
        
        # Clean model matrix
        emptyvars <- names(which(apply(Z_PC, 2, var) == 0))
        Z_PC_clean <- Z_PC[, !colnames(Z_PC) %in% emptyvars]
        
        # Extract PCs
        pcaOut <- prcomp(Z_PC_clean, scale = TRUE, retx = TRUE)
        
        ## Compute and store the cumulative proportion of variance explained by
        ## the component scores:
        rSquared <- cumsum(pcaOut$sdev^2) / sum(pcaOut$sdev^2)

        ## Extract the principal component scores:
        Z_pca <- pcaOut$x[, rSquared <= parms$PCA_pcthresh]

      }
      
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
                                 # printFlag = FALSE,
                                 ridge  = cond$ridge,
                                 method = methods)
      end.time <- Sys.time()
      
      # Identification of models
      imp_PCA_diff_comps <- list(nobs = colSums(!is.na(Z_input[, 1:10])), 
                                 npcs = ncol(imp_PCA_mids$predictorMatrix)-1)
      
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
                  dtIN = Z_out,
                  diff_comp = imp_PCA_diff_comps))
      
      ### END TRYCATCH EXPRESSION
    }, error = function(e){
      err <- paste0("Original Error: ", e)
      print(err)
      return(list(dats = NULL,
                  imps = NULL,
                  time = NULL,
                  mids = NULL,
                  dtIN = NULL,
                  diff_comp = NULL)
      )
    }
    )
    
  } else {
    return(list(dats = NULL,
                imps = NULL,
                time = NULL,
                mids = NULL,
                dtIN = NULL,
                diff_comp = NULL))
  }
}
