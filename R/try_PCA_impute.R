### Title:    imputeHD-comp impute w/ PCA auxiliary variables
### Author:   Edoardo Costantini
### Created:  2019-NOV-26
### Modified: 2019-NOV-26
### Notes:    Main implementation is based on Howard et al 2015
###           - Inclusion of auxiliary vairables: for each variable with missing values, a dataset
###             with only the variable and the first n PC components are use (n should be dedided,
###             look into the paper again)
###           - Initialization: auxiliary variables are initialized with one non stocastich regression
###             run. The idea is that we are not worried about standard error at tgis level. We simply 
###             want a full dataset for the PCA computation.

# load packages
library(tidyverse)
library(gimme)    # for expand.grid.unique

# Prep data ---------------------------------------------------------------

# Create using datagen function
source("./dataGen_test.R")
  set.seed(20191120)
dt <- missDataGen(n=100, p=8)
  dt_c <- dt[[1]] # fully observed
  dt_i <- dt[[2]] # with missings
    dim(dt_i)
    mice::md.pattern(dt_i)

# Define variables with missings
K <- ncol(dt_i)-sum(tail(mice::md.pattern(dt_i),1) == 0) # number of variables needing imputation
K_names <- names(which(colSums(apply(dt_i, 2, is.na)) != 0)) # select the names of these k variables


# Imputation (try my self) ------------------------------------------------

  # Data preparetion
  # Obtain all combinationa and polynomial terms for auxiliary variables
  # FIX; for now it only supports suqared terms
    aux_data <- dt_i[, !names(dt_i) %in% K_names]
    ip_terms <- expand.grid.unique(names(aux_data), names(aux_data), include.equals = T)
    for (t in 1:nrow(ip_terms)) {
      selected <- inter_terms[t,]
      terms_prod <- data.frame(aux_data[, selected[1]]*aux_data[, selected[2]])
      colnames(terms_prod) <- paste0(selected[1], "*", selected[2])
      aux_data <- cbind(aux_data, terms_prod)
    }
    
    esp <- 3
    for (v in 1:ncol(aux_data)) {
      v <- 1
      rep(names(aux_data)[v], esp)
    }
    
  # PCs extraction
    pr_out <- prcomp(aux_data, scale = TRUE)
    pr_out$x
    pr_var <- pr_out$sdev**2
    p_pr_var <- pr_var/sum(pr_var) # proportion of variance explained by component
    PCs_index<- round(cumsum(p_pr_var), 1) <= .5 # selects PCs that approximately explain 50% of the variabce together
  
  # Inlcude auxiliary PCs in the imputation
    dt_i_PCaux <- cbind(dt_i[,1:5], pr_out$x[,PCs_index])
    
  # Impute with mice
    imp <- mice::mice(dt_i[,1:ncol(dt_i)], method = "norm.nob", m = 5)
    imputedList <- complete(imp, "all")
    lapply(imputedList, head)

# Imputation w/ PcAux package ---------------------------------------------
  library(PcAux)
  
  ## Exampole ##
  # Get to know the package with a p < n dataset
    data(iris2)
    md.pattern(iris2)
    
    # First, load and prepare your data:
    cleanData <- prepData(rawData   = iris2,
                          moderators = NULL,           # names of any moderator variables to include 
                                                       # in the initial, single imputation model
                                                       # (if interactType = 2L in the createPcAux functions
                                                       # these would also be the interactions included for that PCA)
                          nomVars   = "Species",       # nominal variables (R: unodered factors)
                          ordVars   = "Petal.Width",   # ordinal variables (R: ordered factors)
                          idVars    = "ID",
                          dropVars  = "Junk",          # list variables you want to get rid of
                          groupVars = "Species")
    
    # Next, create a set of principal component auxiliary variables:
    pcAuxOut <- createPcAux(pcAuxData = cleanData,
                            interactType = 1L, # any of the values in 0, 1, 2, 3 
                                               # indicating the interaction order
                                               # (see detials)
                            maxPolyPow = 3L,   # maximum power used when constructing 
                                               # the polynomial terms
                            nComps    = c(3, 2))
    
    # use PC auxiliaries as the predictors in a multiple imputation run:
    miOut <- miWithPcAux(rawData   = iris2,
                         pcAuxData = pcAuxOut,
                         nImps     = 5)
    getImpData(miOut)
    
    # or work directly with the PC auxiliaries by appending them to the original data
    outData <- mergePcAux(pcAuxData = pcAuxOut, rawData = iris2)
    
    # and then only the PC aux will be used for imputation with mice
    predMat <- makePredMatrix(mergedData = outData)
    
  ## With highdimensional data ##
    source("./dataGen_test.R")
      set.seed(20191120)
    dt <- missDataGen(n=50, p=6)
      dt_c <- dt[[1]] # fully observed
      dt_i <- dt[[2]] # with missings
      dim(dt_i)
    mice::md.pattern(dt_i)
    
    # First, load and prepare your data:
    cleanData <- prepData(rawData   = dt_i,
                          moderators = c(names(dt_i))) # all interactions
    
    # Note on moderators: by inlcuding all variable names in the moderator arguments,
    # you obtain the following
    varCombs <- combn(names(dt_i), 2, simplify = FALSE)
    interact <- data.frame(
      lapply(varCombs,
             function(x, dat) dat[ , x[1]] * dat[ , x[2]],
             dat = dt_i)
    )
    ncol(interact)
    
    # Next, create a set of principal component auxiliary variables:
    pcAuxOut <- createPcAux(pcAuxData = cleanData,
                            interactType = 2L,
                            maxPolyPow = 3L,
                            nComps    = c(3, 2))
    # Imputations
    miOut <- miWithPcAux(rawData   = dt_i,
                         pcAuxData = pcAuxOut,
                         nImps     = 5)
    getImpData(miOut)
    head(dt_i, 9)
    
    # or work directly with the PC auxiliaries by appending them to the original data
    outData <- mergePcAux(pcAuxData = pcAuxOut, rawData = dt_i)
  
    # and then only the PC aux will be used for imputation with mice
    predMat <- makePredMatrix(mergedData = outData)
    
    
    