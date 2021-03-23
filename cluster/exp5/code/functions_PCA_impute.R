### Title:    imputeHD: Functions for PCA imputation strategy
### Author:   Edoardo Costantini
### Created:  2020-05-25
### Notes: Rip off PcAux package functions to make some changes

# PCA manual strategy functions -------------------------------------------

missCheck <- function(x) 
{
  if(missing(x)) return(TRUE)
  if(!is.object(x)) {
    if(length(x) == 0) return(TRUE)
    if(length(x) == 1 & is.null(x) || is.na(x) || x == "") return(TRUE)
  }
  FALSE
}

countLevels <- function(x) 
{
  length(unique(na.omit(x)))
}

codeNomVars = function(data, nomVars)
{
  # Internals
  # data = Z
  # nomVars = colnames(data)[1:2]
  
  "Dummy code nominal factors"
  noms <- colnames(data)[colnames(data) %in% nomVars]
  
  ## Store factor representations:
  facNoms           <- data.frame(data[ , noms])
  colnames(facNoms) <- noms
  dumNoms <- NULL
  
  for(v in noms) {
    ## Expand factors into dummy codes:
    
    ## Make sure missing values are retained in dummy codes:
    oldOpt <- options(na.action = "na.pass")
    
    ## Create/store dummy codes:
    tmp <- data.frame(
      model.matrix(as.formula(paste0("~", v)), data = data)
    )
    dumNames <- colnames(tmp)
    
    ## Reset the na.action option:
    options(na.action = oldOpt$na.action)
    
    ## Remove dummy codes for empty factor levels:
    levVec                 <-  sapply(tmp, countLevels)
    dumNoms[[v]]           <- data.frame(tmp[ , levVec > 1])
    colnames(dumNoms[[v]]) <- dumNames[levVec > 1]
  }
  ## Save dummy code names:
  dumVars <- as.character(unlist(lapply(dumNoms, colnames)))
  
  return(list(dumNoms = dumNoms,
              dumVars = dumVars,
              facNoms = facNoms))
  
}

castOrdVars = function(data, ordVars, toNumeric = TRUE)
{
  
  # Internals
  # data = Z
  # ordVars = colnames(data)[1:2]
  # toNumeric = TRUE
  
  "Cast ordinal factors to numeric variables"
  ## Find ordinal variables that are still on the data set:
  ords <- colnames(data)[colnames(data) %in% ordVars]
  
  if(toNumeric) {
    ## Cast the ordinal variables as numeric:
    if(length(ords) > 1)
      data[ , ords] <-
        data.frame(lapply(data[ , ords], as.numeric))
    else
      data[ , ords] <- as.numeric(data[ , ords])
  }
  else {
    ## Cast back to ordered factors:
    if(length(ords) > 1)
      data[ , ords] <-
        data.frame(lapply(data[ , ords], as.ordered))
    else
      data[ , ords] <- as.ordered(data[ , ords])
  }
  
  return(data)
  
}

castNomVars = function(data, nomVars, moderators, toNumeric = TRUE)
{
  
  # Internals
  # data = Z
  # nomVars = colnames(data)[1:2]
  # moderators = colnames(data)[1]
  # toNumeric = TRUE
  
  "Swap factor and dummy-coded representations of nominal variables"
  
  dumObj <- codeNomVars(data = data, nomVars = nomVars)
  dumNoms <- dumObj$dumNoms
  dumVars <- dumObj$dumVars
  facNoms <- dumObj$facNoms
    
  if(toNumeric) {# Replace factors in the data with dummy codes
    otherNames     <-  setdiff(colnames(data), nomVars)
    
    
    data           <- data.frame(data[ , otherNames],
                                  do.call(cbind, dumNoms)
    )
    colnames(data) <- c(otherNames, dumVars)
    
    ## Update the moderators with dummy code names:
    check <- moderators %in% nomVars
    if(any(check)) {
      frozenMods <- moderators
      moderators <-
        c(moderators[!check],
          unlist(lapply(dumNoms[moderators[check]], colnames),
                 use.names = FALSE)
        )
    }
  }
  else {# Undo dummy coding
    data <-
      data.frame(data[ , setdiff(colnames(data), dumVars)],
                 facNoms)
    
    ## Revert to original moderator list:
    moderators <- frozenMods
  }
  
  return(data)
  
}

computeInteract = function(data, 
                           idVars, ordVars, nomVars,
                           moderators,
                           intMeth = 1L)
{
  
  # Internals
  # data = Z_aux
  # idVars = colnames(data)
  # ordVars = colnames(data)[1:3]
  # nomVars = colnames(data)[4:5]
  # moderators = colnames(data) #colnames(data)[8:10]
  # intMeth = 1L
  
  "Calculate interaction terms"
  # if(length(lin) > 0) # Do we have linear PcAux?
  #   pcNames <- setdiff(colnames(lin), idVars)
  
  ## Cast factors to numeric:
  if(length(ordVars) > 0) 
    data <- castOrdVars(data = data, ordVars = ordVars)
  
  if(length(nomVars) > 0) 
    data <- castNomVars(data = data, 
                        moderators = moderators,
                        nomVars = nomVars)
  
  ## Define moderators and focal predictors:
  if(intMeth < 3) mods <- moderators
  else            mods <- colnames(data)
  
  if(intMeth == 1) focal <- setdiff(colnames(data), mods)
  else             focal <- pcNames
  
  ## Generate a list of interacted variable pairs:
  if(length(focal) == 0) {
    ## All possible two-way interactions:
    varCombs <- combn(mods, 2, simplify = FALSE)
  }
  else {
    ## Subset of interactions defined by user:
    varCombs <- list()
    for(m in mods)
      for(f in focal)
        varCombs[[paste0(m, f)]] <- c(m, f)
  }
  
  ## Generate variable names for interaction terms:
  intVars <- as.character(sapply(varCombs, paste0, collapse = "."))
  
  ## Compute the interaction terms:
  if(intMeth == 1)
    interact <- data.frame(
      lapply(varCombs,
             function(x, dat) dat[ , x[1]] * dat[ , x[2]],
             dat = data)
    )
  else
    interact <<- data.frame(
      lapply(varCombs,
             function(x, dat, pc) dat[ , x[1]] * pc[ , x[2]],
             dat = data,
             pc  = lin[ , pcNames])
    )
  colnames(interact) <- intVars
  
  ## Remove dummy codes for empty cells:
  levVec   <- sapply(interact, countLevels)
  interact <- interact[ , levVec > 1]
  intVars  <- colnames(interact)
  
  ## Cast dummy codes as factors:
  dumFlag <-
    sapply(interact,
           function(x) all(unique(na.omit(x)) %in% c(0, 1))
    )
  
  if(sum(dumFlag) > 1) {
    interact[ , dumFlag] <-
      data.frame(lapply(interact[ , dumFlag], as.factor))
  } else {
    if(sum(dumFlag) == 1)
      interact[ , dumFlag] <- as.factor(interact[ , dumFlag])
  }
  
  ## If any PcAux are involved, orthogonalize the interaction terms
  ## w.r.t. the linear PcAux scores:
  # if(intMeth > 1)        
  #   for(v in 1 : ncol(interact))
  #     interact[ , v] <<-
  #   .lm.fit(y = interact[ , v],
  #           x = as.matrix(lin[ , pcNames]))$resid
  
  return(interact)
}

computePoly = function(data, ordVars, nomVars, maxPower = 2L)
{
  # For internals
  # data <- Z_aux
  # nomVars <- colnames(data)[1:3] # will be excluded from polyterms
  # ordVars <- colnames(data)[5:10] # will be excluded from polyterms
  # maxPower <- 3
  
  "Compute polynomial terms"
  ## Cast ordered factors to numeric:
  if(length(ordVars) > 0) 
    data <- castOrdVars(data = data, ordVars = ordVars)
  
  dataNames <- setdiff(colnames(data), nomVars)
  
  ## Construct a formula to define the polynomial transformations:
  form <- as.formula(
    paste0("~",
           paste0("I(",
                  paste(dataNames,
                        rep(c(2 : maxPower),
                            each = length(dataNames)),
                        sep = "^"),
                  ")",
                  collapse = " + "
           )
    )
  )
  
  ## Make sure missing values are retained in dummy codes:
  oldOpt <- options(na.action = "na.pass")
  
  ## Create the polynominal terms:
  if(length(dataNames) == 1) # Hack for only one target variable
    poly <<- data.frame(
      model.matrix(form,
                   data = as.data.frame(          
                     list(data[ , dataNames]),
                     col.names = dataNames)       
      )[ , -1]
    )
  else
    poly <- data.frame(
      model.matrix(form, data = data[ , dataNames])
    )[ , -1]
  
  ## Reset the na.action option:
  options(na.action = oldOpt$na.action)
  
  return(poly)
  
}



# Junk --------------------------------------------------------------------

doPCA <- function(lin) {
  # For internals
  lin = 0
  nComps = c(5, 0)
  intMeth = 1L
  maxPower = 2L
  
  ## Are we extracting linear or nonlinear PC scores?
  if(length(map$pcAux$lin) == 0) {linVal <- "lin";    pcType <- 1}
  else                           {linVal <- "nonLin"; pcType <- 2}
  
  ## Do we need to parse the nComps argument?
  parseCheck <- is.infinite(map$nComps[pcType]) |
    (map$nComps[pcType] < 1 & map$nComps[pcType] != 0)
  
  if(linVal == "lin") {
    # if(map$verbose > 0) # alqways verbose
    cat("\nCalculating linear principal component scores...\n")
    
    if(map$intMeth != 1 & map$maxPower > 1) {
      ## Replace the raw polinomial values with the imputed versions from
      ## the data:
      map$poly <- map$data[ , colnames(map$poly)]
      ## Remove polynomial terms from the data before extracting linear
      ## PcAux:
      map$data <- with(map,
                       data[ , setdiff(colnames(data), colnames(poly))]
      )
    }
    
    ## Cast ordinal factors to numeric formats:
    if(!missCheck(map$ordVars)) map$castOrdVars() 
    ## Dummy code nominal variables:
    if(!missCheck(map$nomVars)) map$castNomVars() 
    
    # ## Cast any remaining factors to numeric formats:
    # check <- sapply(map$data, is.factor)
    # if(sum(check) > 1)
    #   map$data[ , check] <- data.frame(lapply(map$data[ , check], f2n))
    # if(sum(check) == 1)
    #   map$data[ , check] <- f2n(map$data[ , check])
    # 
    # if(!map$simMode & !parseCheck) {
    #   ## Make sure the number of PC scores we want is less than the number
    #   ## of columns in our data object:
    #   if(map$nComps[1] > ncol(map$data)) {
    #     warnFun("linPcNum", map)
    #     map$nComps[1] <- ncol(map$data)
    #   }
    # }
  } else {# We already have linear component scores
    # if(map$verbose > 0) # always verbose
    cat("\nCalculating nonlinear principal component scores...\n")
    
    ## Redefine the data object:
    if(map$intMeth > 1) map$data <- map$interact
    
    if(map$maxPower > 1) {
      ## Orthogonalize the polynomial terms w.r.t. the linear PcAux:
      map$poly <- data.frame(
        lapply(map$poly,
               function(y, pc)
                 .lm.fit(y = y, x = as.matrix(pc))$residuals,
               pc = with(map,
                         pcAux$lin[ , setdiff(colnames(pcAux$lin),
                                              idVars)]
               )
        )
      )
      if(map$intMeth > 1) map$data <- data.frame(map$data, map$poly)
      else                map$data <- map$poly
    }
    colnames(map$data) <- c(colnames(map$interact), colnames(map$poly))
    
    ## Remove the contents of the 'interact' and 'poly' fields when they are
    ## no longer necessary:
    map$interact <- "Removed to save resources"
    map$poly     <- "Removed to save resources"
    
    if(!map$simMode & !parseCheck) {
      ## Make sure the number of PC scores we want is less than
      ## the number of columns in our data object:
      if(map$nComps[2] > ncol(map$data)) {
        warnFun("nonLinPcNum", map)
        map$nComps[2] <- ncol(map$data)
      }
    }
  }# CLOSE if(length(map$pcAux$lin) == 0)
  
  ## Execute the principal component analysis:
  if(map$pcaMemLev == 0) {
    ## Higher numerical accuracy, but more memory usage
    pcaOut <- prcomp(map$data, scale = TRUE, retx  = TRUE)
    
    ## Compute and store the cumulative proportion of variance explained by
    ## the component scores:
    map$rSquared[[linVal]] <- cumsum(pcaOut$sdev^2) / sum(pcaOut$sdev^2)
    
    ## Set component counts when some are defined in terms of variance
    ## explained:
    if(parseCheck) map$setNComps(type = pcType)
    
    ## Extract the principal component scores:
    if(length(map$idCols) == 0)
      map$pcAux[[linVal]] <- pcaOut$x[ , 1 : map$nComps[pcType]]
    else
      map$pcAux[[linVal]] <-
      data.frame(map$idCols, pcaOut$x[ , 1 : map$nComps[pcType]])
  } else if(map$pcaMemLev == 1) {
    ## Save memory at the expense of numerical accuracy
    pcaOut <-
      simplePca(map = map, lv = pcType, parse = parseCheck, scale = TRUE)
  } else {
    errFun("badPcaMemLev", map = map)
  }
  
  ## Remove the contents of the 'data' field when they're no longer needed:
  if(linVal == "nonLin") map$data <- "Removed to save resources"
  
  ## Give some informative column names:
  colnames(map$pcAux[[linVal]]) <-
    c(map$idVars,
      paste0(ifelse(linVal == "lin", "linPC", "nonLinPC"),
             c(1 : map$nComps[ifelse(linVal == "lin", 1, 2)])
      )
    )    
  if(map$verbose > 0) cat("Complete.\n")
}
