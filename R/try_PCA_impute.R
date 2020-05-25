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

rm(list = ls())

# load packages
library(tidyverse)
library(gimme)    # for expand.grid.unique
library(PcAux)
library(mice)

# Prep data ---------------------------------------------------------------

# Create using datagen function
source("./dataGen_test.R")
  set.seed(20191120)
dt <- missDataGen(n=100, p=10)
  dt_c <- dt[[1]] # fully observed
  Z <- dt[[2]][,-c(1:3)] # with missings
    dim(Z)
    mice::md.pattern(Z)

# Define variables with missings
K <- ncol(Z)-sum(tail(mice::md.pattern(Z),1) == 0) # number of variables needing imputation
K_names <- names(which(colSums(apply(Z, 2, is.na)) != 0)) # select the names of these k variables


# Imputation (try my self) ------------------------------------------------
  ycont_formula <- c("yI ~ V1 + V2 + V3 + V4") # some model of interest
  # Select auxiliary variables (PCA target)
  Z_aux <- Z[, -which(sapply(names(Z), grepl, x = ycont_formula))]
  
  # 1. Single Imputation of auxiliary vairbales if missing
  Z_aux_mids <- mice(Z_aux, m = 1, maxit = 1, method = "norm.nob")
  Z_aux_fil <- complete(Z_aux_mids)
  
  # 2. Create Set of Auxiliary vairables interactions and polynomials
  varCombs <- combn(names(Z_aux), 2, simplify = FALSE)
  interact <- data.frame(
    sapply(1:length(varCombs),
           function(x) Z_aux[ , varCombs[[x]][1]] * Z_aux[ , varCombs[[x]][2]])
  )
  
  polynom2 <- sapply(Z_aux, function(x) x^2)
  
  Z_aux <- cbind(Z_aux, interact, polynom2)
  
  # Extract PCs
  ncol(Z_aux)
  pr.out <- prcomp(Z_aux, scale = TRUE)
  exp_var <- pr.out$sdev^2/sum(pr.out$sdev^2)*100
  Z_pca <- pr.out$x[, cumsum(exp_var) < 50]
  Z_mod <- Z[, which(sapply(names(Z), grepl, x = ycont_formula))]
  
  # Impute
  imp_PCA_mids <- mice::mice(cbind(Z_mod, Z_pca),
                             m = 5,
                             maxit = 10,
                             ridge = 1e-5,
                             method = "norm")
  imputedList <- complete(imp, "all")
  lapply(imputedList, head)

# Imputation w/ PcAux package ---------------------------------------------
  library(PcAux)
  
  # Prepare all-variables data set
  # The step (prepData) involves:
  # 1: castData: cast variables to declared types
  #             a) Flag variable types
  #             b) Cast all variables to the appropriate measurement level
  #             c) Dummy code nominal variables if in creatPcAux()
  #             d) Center continuous variables (not scaled)
  # 2. cleanData: Find and (possibly) remove problematic data columns 
  #               (i.e., variables with few or no observations/constants)
  # 3. findCollin: Flag variables with perfect bivariate correlations
  # Step 2 and 3 are skipped in simulation mode (simMode = TRUE)
  cleanData <- prepData(rawData   = Z,
                        moderators = NULL, 
                         # moderators for initial single imputation
                         # NULL makes it everything
                        simMode = TRUE, # this measn that only cast data is applyed
                        verbose = TRUE) 
  cleanData$data
  
  # Create PCs
  # 1. pcAuxData$computeInteract(): Compute interactions for use during initial imputation:
  # 2. pcAuxData$computePoly(): Compute polynomials for use during initial imputation
  # 3. doSingleImputation(map = pcAuxData): Execute the initial, single imputation
  # 4. doPCA(map = pcAuxData): Extract the linear principal component scores
  #     a) extracts linear components (i.e. polynomial terms are exluded before extraction)
  #         - ordinal factors are cast to numeric
  #         - nominal are cast to dummies
  #         - prcomp(targetData, scale = TRUE, retx  = TRUE)
  #         - extract component scores based on variance or number
  #     b) extracts non-linear components if required
  
  pcAuxOut <- createPcAux(pcAuxData = cleanData,
                          interactType = 1L, 
                            # inclusion of itneractions int eh intial single imputation
                            # all two-way interactions between the observed variables 
                            # and the variables specified in the moderators argument 
                            # of prepData
                          maxPolyPow = 3L,   # maximum power used when constructing 
                          # the polynomial terms
                          nComps    = c(5, 5))
  
  # Perform imputations of interest
  # 1. mergePcAux: Combine principal component auxiliaries with raw data
  PCaux_Data <- mergePcAux(pcAuxData = pcAuxOut,
                           rawData   = Z,
                           nComps    = c(5, 0))
  # 2. construct predictor matrix for mice: 
    nLin    <- length(grep("^linPC\\d", colnames(PCaux_Data)))
    nNonLin <- length(grep("^nonLinPC\\d", colnames(PCaux_Data)))
    
    predMat <- makePredMatrix(mergedData   = PCaux_Data,
                              nLinear      = nLin,
                              nNonLinear   = nNonLin)
  # 3. Impute data
    mice(PCaux_Data,
         m               = 5,
         maxit           = 1L,
         predictorMatrix = predMat)
    
  # Use all of it in a function
  miOut <- miWithPcAux(rawData   = Z,
                       pcAuxData = pcAuxOut, # object produced by createPcAux
                       nImps     = 5,
                       nComps = c(5, 0))
  

# Real data examples ------------------------------------------------------
# > Iris data ####
  data(iris2, package = "PcAux")
  dim(iris2)
  md.pattern(iris2, plot = FALSE)
  
  # First, load and prepare your data:
  cleanData <- prepData(rawData   = iris2,
                        moderators = NULL,           # names of any moderator variables to include 
                                                     # in the initial, single imputation model
                                                     # (if interactType = 2L in the createPcAux functions
                                                     # these would also be the interactions included for that PCA)
                                                     # NULL includes every variable as moderator
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
                          nComps    = c(.6, .6))
  
  PCaux_Data <- mergePcAux(pcAuxData = pcAuxOut,
                           rawData   = iris2,
                           nComps    = c(.6, .6))
  
  pcAuxOut
  str(pcAuxOut)
  pcAuxOut$pcAux
  
  # use PC auxiliaries as the predictors in a multiple imputation run:
  miOut <- miWithPcAux(rawData   = iris2,
                       pcAuxData = pcAuxOut,
                       nImps     = 5)
  getImpData(miOut)
  miOut$miceObject$predictorMatrix # only the pcs are used to impute
  # or work directly with the PC auxiliaries by appending them to the original data
  outData <- mergePcAux(pcAuxData = pcAuxOut, rawData = iris2)
  
  # and then only the PC aux will be used for imputation with mice
  predMat <- makePredMatrix(mergedData = outData)
