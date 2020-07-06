### Title:    imputeHD-comp impute with sequential CART
### Author:   Edoardo Costantini
### Created:  2019-DEC-10
### Modified: 2019-DEC-16
### Notes:    Main implementation based on Shah et al. This script serves the sole purpose
###           of helping me understand what goes on in the packages CALIBERrfimpute and mice when
###           imputing with random forests. It is not meant to be used for other purposes.

# load packages
library(tidyverse)
library(mice)
library(randomForest)
library(CALIBERrfimpute) # package from Shah et al 2014

# Load all purpose functions
  source("./functions_allpurp.R")

# data --------------------------------------------------------------------
  # Create using datagen function
  source("./dataGen_test.R")
  set.seed(20191213)
  dt <- missDataGen(n=1e3, p=20)
    dt_c <- dt[[1]] # fully observed
    dt_i <- dt[[2]] # with missings
    dim(dt_i)
  
# Functions study ---------------------------------------------------------
  # Prep data
    K <- ncol(dt_i)-sum(tail(mice::md.pattern(dt_i),1) == 0) # number of variables needing imputation
    K_names <- names(which(colSums(apply(dt_i, 2, is.na)) != 0)) # select the names of these k variables
    Z <- dt_i    # dataset with missing values
    p <- ncol(Z) # number of variables [INDEX with j]
    n <- nrow(Z) # number of observations
    colnames(Z) <- paste0(rep("z_", p), seq(1, p))
    r <- colSums(apply(Z, 2, function (x) {!is.na(x)} )) # vector with number of observed values in a j variable

    # Define variables with missings
    l <- ncol(Z)-sum(tail(mice::md.pattern(Z),1) == 0) # number of variables needing imputation
    l_names <- names(which(colSums(apply(Z, 2, is.na)) != 0)) # select the names of these k variables
    
    # Define a dataset to recieve initialized missing values
    Z0 <- Z 
    
    # Define names of variables w/ missing values (general and by measurment scale)
    vartypes <- rbind(lapply(lapply(Z0[, names(Z0) %in% l_names], class), paste, collapse = " "))
    contVars <- colnames(vartypes)[vartypes == "numeric"] # names of continuous variables for selection
    factVars <- colnames(vartypes)[vartypes == "factor"] # includes factors and ordered factors
    ordeVars <- colnames(vartypes)[vartypes == "ordered factor"] # includes factors and ordered factors
    
    # Make oredered factors as numeric
    Z0[, ordeVars] <- as.numeric(Z0[, ordeVars])
    
    # Impute sample means for continuous variables and odered factors
    s_means <- apply(Z0[, c(contVars, ordeVars)], 2, mean, na.rm = TRUE) # sample means
    for (j in 1:length(c(contVars, ordeVars))) {
      Z0 <- Z0 %>% mutate_at(vars(c(contVars, ordeVars)[j]),
                             ~replace(., is.na(.), s_means[j])
      )
    }
    # Impute most common level for unordered factors
    for (j in 1:length(factVars)) {
      x <- addNA(Z0[, factVars[j]])
      m_commo <- names(which.max(table(x)))
      levels(x) <- c(levels(Z0[, factVars[j]]), 99)
      x[x == 99] <- m_commo
      x <- droplevels(x)
      Z0[, factVars[j]] <- x
    }
    
    Zm <- Z0
    
    
#### > mice.impute.rfcont() ####
    # Which variable to impute
    k <- 4 # for imputing one variable ()
    
    # inputs
    y <- Z[, l_names[k]]
    ry <- !is.na(y) # response variable T = observed, F = NA
    x <- Zm[, -which(colnames(Zm) == l_names[k])]
    ntree_cont    = 7    # number of trees
    nodesize_cont = NULL # minimum size of nodes
    maxnodes_cont = NULL # maximum number of nodes
    ntree = NULL         # 
    
    # function
    #x <- as.matrix(x) # mixed data does not become matrix easily
    bootsample <- sample(x=sum(ry), replace = TRUE) #sample takes place from 1:x
    
    # Random forest parameters
    if (is.null(ntree_cont) & !is.null(ntree)) {
      ntree_cont <- ntree
    }
    if (is.null(ntree_cont)) {
      if (is.null(getOption("CALIBERrfimpute_ntree_cont"))) {
        ntree_cont <- 10
      }
      else {
        ntree_cont <- getOption("CALIBERrfimpute_ntree_cont")
      }
    }
    if (is.null(nodesize_cont)) {
      if (is.null(getOption("CALIBERrfimpute_nodesize_cont"))) {
        nodesize_cont <- 5
      }
      else {
        nodesize_cont <- getOption("CALIBERrfimpute_nodesize_cont")
      }
    }
    if (is.null(maxnodes_cont)) {
      maxnodes_cont <- getOption("CALIBERrfimpute_maxnodes_cont")
    }
    if (ntree_cont > 1) { # SELECTION of BOOTSTRAPPED OBSERVATIONS
      yobs <- y[ry][bootsample]
      xobs <- x[ry, , drop = FALSE][bootsample, , drop = FALSE]
    } else {
      yobs <- y[ry]
      xobs <- x[ry, , drop = FALSE]
    }
    
    # Actual imputation
    xmiss <- x[!ry, , drop = FALSE]
    rf <- randomForest(xobs, yobs, ntree = ntree_cont, nodesize = nodesize_cont, 
                       maxnodes = maxnodes_cont)
    
    yhat <- predict(rf, xmiss,
                    predict.all=FALSE) # as means of the conditonal distributions
                                      # the values predicted using all of the random forest
                                      # trees are used (aggregate). If predict.all was FALSE
                                      # then, both aggregate and individual predictions (done
                                      # by each of the trees) would be reported.
    
    sqrtOOB <- sqrt(rf$mse[ntree_cont]) # for the last tree (just one of the 
                                        # ones we could have selected) the mean squared error of 
                                        # prediction made on the train data xobs
                                        # This does not make sense to me! Why the last tree? <- ????????????
    yimp <- rnorm(length(yhat), mean = yhat, sd = sqrtOOB)
    #return(yimp)
    
#### > mice.impute.norm.boot ####
    # Shah et al implementation of RF is said to be same as mice.impute.norm.boot
    # but with the further bootstrap sampling implied by the RF (of the variables to be used)
    function (y, ry, x, wy = NULL, ...) 
    {
      if (is.null(wy)) 
        wy <- !ry
      x <- cbind(1, as.matrix(x))
      n1 <- sum(ry)
      s <- sample(n1, n1, replace = TRUE)
      ss <- s
      dotxobs <- x[ry, , drop = FALSE][s, ]
      dotyobs <- y[ry][s]
      p <- estimice(dotxobs, dotyobs, ...)
      sigma <- sqrt((sum(p$r^2))/(n1 - ncol(x) - 1))
      return(x[wy, ] %*% p$c + rnorm(sum(wy)) * sigma)
    }
      
#### > mice.impute.rfcat() ####
    # Which variable to impute
    k <- 1 # for imputing one variable ()
    
    # inputs
    y <- Z[, l_names[k]]
    ry <- !is.na(y) # response variable T = observed, F = NA
    x <- Zm[, -which(colnames(Zm) == l_names[k])]
    ntree_cat = NULL
    nodesize_cat = NULL
    maxnodes_cat = NULL
    ntree = NULL
    
    if (is.logical(y)) {
      convertlogical <- TRUE
      y <- as.factor(y)
    } else {
      convertlogical <- FALSE
    }
    # Random forest parameters
    if (is.null(ntree_cat) & !is.null(ntree)) {
      ntree_cat <- ntree
    }
    if (is.null(ntree_cat)) {
      if (is.null(getOption("CALIBERrfimpute_ntree_cat"))) {
        ntree_cat <- 10
      }
      else {
        ntree_cat <- getOption("CALIBERrfimpute_ntree_cat")
      }
    }
    if (is.null(nodesize_cat)) {
      if (is.null(getOption("CALIBERrfimpute_nodesize_cat"))) {
        nodesize_cat <- 1
      }
      else {
        nodesize_cat <- getOption("CALIBERrfimpute_nodesize_cat")
      }
    }
    if (is.null(maxnodes_cat)) {
      maxnodes_cat <- getOption("CALIBERrfimpute_maxnodes_cat")
    }
    
    # Bootstrapping
    bootsample <- sample(sum(ry), replace = TRUE)
    yobs <- y[ry][bootsample]
    xobs <- x[ry, , drop = FALSE][bootsample, , drop = FALSE]
    xmiss <- x[!ry, , drop = FALSE]
    
    # Some stuff to fix the labeling (not imp)
    missinglevels <- (table(yobs) == 0)
    newlevels <- rep(NA_integer_, length(levels(y)))
    newlevels[!missinglevels] <- 1:sum(!missinglevels)
    labels <- levels(y)
    oldlevels <- 1:length(levels(y))
    changes <- !identical(newlevels, 1:length(levels(y)))
    if (changes) {
      temp <- data.frame(id_yobs = 1:length(yobs), fac = as.integer(yobs))
      lookup <- data.frame(fac = oldlevels, new = factor(newlevels))
      temp <- merge(temp, lookup, all.x = TRUE)
      yobs <- temp[order(temp$id_yobs), "new"]
    }
    
    # perform the imputations
    trees <- lapply(1:ntree_cat, function(x) {
      if (length(unique(yobs)) == 1) {
        yobs[1]
      }
      else {
        randomForest(xobs, yobs, ntree = 1, nodesize = nodesize_cat, 
                     maxnodes = maxnodes_cat)
      }
    })
    
    # This does not work beacuse you are keeping the data.frame structure instead
    # of the matrix form.
    yimp <- apply(xmiss, 1, function(x) {
      thetree <- trees[[sample(ntree_cat, 1)]]
      if ("randomForest" %in% class(thetree)) {
        predict(thetree, x)
      }
      else {
        thetree
      }
    })
    # This is what should happen inside for all the rows of xmiss
    thetree <- trees[[sample(ntree_cat, 1)]] # random tree selected and used for prediction
    predict(thetree, xmiss[1,])
    
    if (changes) {
      temp <- data.frame(id_yimp = 1:length(yimp), fac = as.integer(yimp))
      lookup <- data.frame(fac = newlevels, old = factor(oldlevels))
      levels(lookup$old) <- labels
      temp <- merge(temp, lookup, all.x = TRUE)
      yimp <- temp[order(temp$id_yimp), "old"]
    }
    if (convertlogical) {
      yimp <- as.logical(yimp == "TRUE")
    }
    return(yimp)
    

    
#### > mice() ####
    files.sources = list.files("../notes/code_notes/mice-master/R")
    sapply(paste0("../notes/code_notes/mice-master/R/", files.sources), source)
    # load all functionms from the mice package (for helper functions like the checks)
    
# force the imputation methods of CALIBER package in place of methods 
# such as mice.impute.norm.boot or mice.impute.rf

      set.seed(20191216)
      data = dt_i
      m = 5
      method = "sh.rf"
      predictorMatrix = make.predictorMatrix(data)
      blocks = NULL
      where = NULL
      visitSequence = NULL
      post=NULL
      blots=NULL
      data.init = NULL
      maxit = 5
      printFlag = TRUE

# mice <- function(data, m = 5, 
#                  method = NULL,
#                  predictorMatrix,
#                  where = NULL,
#                  blocks,
#                  visitSequence = NULL,
#                  formulas,
#                  blots = NULL,
#                  post = NULL,
#                  defaultMethod = c("pmm", "logreg", "polyreg", "polr"),
#                  maxit = 5, printFlag = TRUE, seed = NA,
#                  data.init = NULL,
#                  ...) {
# call <- match.call()
# check.deprecated(...)
# if (!is.na(seed)) set.seed(seed)
  
  # check form of data and m
  data <- check.dataform(data)
  m <- check.m(m)
  
  # determine input combination: predictorMatrix, blocks, formulas
  mp <- missing(predictorMatrix)
  mb <- FALSE
  # mf <- missing(formulas) # alternative to predictorMatrix
  mf <- FALSE
  
  # # case A
  # if (mp & mb & mf) {
  #   # blocks lead
  #   blocks <- make.blocks(colnames(data))
  #   predictorMatrix <- make.predictorMatrix(data, blocks)
  #   formulas <- make.formulas(data, blocks)
  # }
  # # case B
  # if (!mp & mb & mf) {
    # predictorMatrix leads (assuming this is the case, nothing else is important)
    predictorMatrix <- check.predictorMatrix(predictorMatrix, data)
    blocks <- make.blocks(colnames(predictorMatrix), partition = "scatter")
    formulas <- make.formulas(data, blocks, predictorMatrix = predictorMatrix)
  # }
  # 
  # # case C
  # if (mp & !mb & mf) {
  #   # blocks leads
  #   blocks <- check.blocks(blocks, data)
  #   predictorMatrix <- make.predictorMatrix(data, blocks)
  #   formulas <- make.formulas(data, blocks)
  # }
  # 
  # # case D
  # if (mp & mb & !mf) {
  #   # formulas leads
  #   formulas <- check.formulas(formulas, data)
  #   blocks <- construct.blocks(formulas)
  #   predictorMatrix <- make.predictorMatrix(data, blocks)
  # }
  # 
  # # case E
  # if (!mp & !mb & mf) {
  #   # predictor leads
  #   blocks <- check.blocks(blocks, data)
  #   z <- check.predictorMatrix(predictorMatrix, data, blocks)
  #   predictorMatrix <- z$predictorMatrix
  #   blocks <- z$blocks
  #   formulas <- make.formulas(data, blocks, predictorMatrix = predictorMatrix)
  # }
  # 
  # # case F
  # if (!mp & mb & !mf) {
  #   # formulas lead
  #   formulas <- check.formulas(formulas, data)
  #   predictorMatrix <- check.predictorMatrix(predictorMatrix, data)
  #   blocks <- construct.blocks(formulas, predictorMatrix)
  # }
  # 
  # # case G
  # if (mp & !mb & !mf) {
  #   # blocks lead
  #   blocks <- check.blocks(blocks, data, calltype = "formula")
  #   formulas <- check.formulas(formulas, blocks)
  #   predictorMatrix <- make.predictorMatrix(data, blocks)
  # }
  # 
  # # case H
  # if (!mp & !mb & !mf) {
  #   # blocks lead
  #   blocks <- check.blocks(blocks, data)
  #   formulas <- check.formulas(formulas, data)
  #   predictorMatrix <- check.predictorMatrix(predictorMatrix, data, blocks)
  # }
  
  chk   <- check.cluster(data, predictorMatrix)  
  where <- check.where(where, data, blocks)
  visitSequence <- check.visitSequence(visitSequence, data = data, 
                                       where = where, blocks = blocks)
  method <- check.method(method = method, data = data, where = where, 
                         blocks = blocks, defaultMethod = defaultMethod)
  post <- check.post(post, data)
  blots <- check.blots(blots, data, blocks)
  
  # data frame for storing the event log
  state <- list(it = 0, im = 0, dep = "", meth = "", log = FALSE)
  loggedEvents <- data.frame(it = 0, im = 0, dep = "", meth = "", out = "")
  
  # edit imputation setup
  setup <- list(method = method,
                predictorMatrix = predictorMatrix,
                visitSequence = visitSequence, 
                post = post)
  setup <- edit.setup(data, setup)
  method <- setup$method
  predictorMatrix <- setup$predictorMatrix
  visitSequence <- setup$visitSequence
  post <- setup$post
  
  # initialize imputations
  nmis <- apply(is.na(data), 2, sum)
  imp <- initialize.imp(data, m, where, blocks, visitSequence, 
                        method, nmis, data.init)
  
  # and iterate...
  from <- 1
  to <- from + maxit - 1
  q <- sampler(data, m, where, imp, blocks, method, visitSequence,  # the sampler function is the one that calls the method
               predictorMatrix, formulas, blots, post, c(from, to), 
               printFlag)
  
  if (!state$log) loggedEvents <- NULL
  if (state$log) row.names(loggedEvents) <- seq_len(nrow(loggedEvents))
  
  ## save, and return
  midsobj <- list(data = data, imp = q$imp, m = m,
                  where = where, blocks = blocks, 
                  call = call, nmis = nmis, 
                  method = method,
                  predictorMatrix = predictorMatrix,
                  visitSequence = visitSequence,
                  formulas = formulas, post = post, 
                  blots = blots,
                  seed = seed, 
                  iteration = q$iteration,
                  lastSeedValue = .Random.seed, 
                  chainMean = q$chainMean,
                  chainVar = q$chainVar, 
                  loggedEvents = loggedEvents,
                  version = packageVersion("mice"),
                  date = Sys.Date())
  oldClass(midsobj) <- "mids"
  
  if (!is.null(midsobj$loggedEvents)) 
    warning("Number of logged events: ", nrow(midsobj$loggedEvents),
            call. = FALSE)
  return(midsobj)
# }