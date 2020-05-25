### Title:    Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-05-19

impute_PCA <- function(Z, parms = parms){
  
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
  # Z = Xy_mis
  
  ## body:
  if(parms$meth_sel$MI_PCA == TRUE){
    start.time <- Sys.time()
    
    # Data
    p <- ncol(Z) # number of variables [INDEX with j]
    n <- nrow(Z) # number of observations
    
    # Separate analysis model vairables from auxiliary variables 
    Z_aux <- Z[, -which(sapply(names(Z), grepl, x = parms$formula))]
    Z_mod <- Z[, which(sapply(names(Z), grepl, x = parms$formula))]
    
    # Variables descriptions
    vartypes <- apply(Z_aux, 2, class)
      contVars <- names(vartypes)[vartypes == "numeric"]
      factVars <- names(vartypes)[vartypes == "factor"]
      ordeVars <- names(vartypes)[vartypes == "ordered factor"]
      
    
    # 1. Single Imputation of auxiliary vairbales if missing 
    # (fully observed axuliaries at the moment)
    # Z_aux_mids <- mice(Z_aux, m = 1, maxit = 1, method = "norm.nob")
    # Z_aux <- complete(Z_aux_mids)
    
    # 2. Create Set of Auxiliary variables interactions and polynomials
    interact <- computeInteract(Z_aux, 
                                idVars = colnames(Z_aux),
                                ordVars = ordeVars,
                                nomVars = factVars,
                                moderators = colnames(Z_aux) )
    
    polynom <- computePoly(Z_aux, 
                           ordVars = ordeVars,
                           nomVars = factVars,
                           maxPower = 2L)
    
    Z_aux <- cbind(Z_aux, interact, polynom)
    
    # Extract PCs
    pcaOut <- prcomp(Z_aux, scale = TRUE, retx  = TRUE)
    
    ## Compute and store the cumulative proportion of variance explained by
    ## the component scores:
    rSquared <- cumsum(pcaOut$sdev^2) / sum(pcaOut$sdev^2)
    
    ## Extract the principal component scores:
    Z_pca   <- pcaOut$x[, rSquared < .75]
    
    # Impute
    imp_PCA_mids <- mice::mice(cbind(Z_mod, Z_pca),
                               m = parms$ndt,
                               maxit = parms$iters,
                               ridge = 1e-5,
                               method = "norm")
    end.time <- Sys.time()
    
    # Store results
    imp_PCA_dats <- mice::complete(imp_PCA_mids, "all")
    imp_PCA_imps <- imp_PCA_mids$imp
    imp_PCA_time <- difftime(end.time, start.time, units = "mins")
    
    return(list(dats = imp_PCA_dats,
                imps = imp_PCA_imps,
                time = imp_PCA_time,
                mids = imp_PCA_mids))
  } else {
    return(list(dats = NULL,
                imps = NULL,
                time = NULL,
                mids = NULL))
  }
}
