### Title:    understanding the mice code
### Author:   Edoardo Costantini
### Created:  2019-DEC-16
### Modified: 2019-DEC-17
### Notes:    This script is copy of the mice function from van Buuren that I used to
###           take notes on what happens inside the function. Not contribution or improvement
###           claimed.

library(mice)

# Load helper functions
files.sources = list.files("../notes/code_notes/mice-master/R")
sapply(paste0("../notes/code_notes/mice-master/R/", files.sources), source)
  # load all functionms from the mice package (for helper functions like the checks)

# Data
source("./dataGen_test.R")
  set.seed(20191213)
  dt <- missDataGen(n=20, p=10)
  dt_i <- dt[[2]] # with missings

## MICE call
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
  
  # Define inputs
  data = dt_i
  m = 5
  method = "rf"
  predictorMatrix = make.predictorMatrix(data)
  blocks = NULL
  where = NULL
  visitSequence = NULL
  post=NULL
  blots=NULL
  data.init = NULL
  maxit = 5
  printFlag = TRUE
  
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
  #
  # ...
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
  
  #### in the function:
  # q <- sampler(data, m, where, imp, blocks, method, visitSequence,  # the sampler function is the one that calls the method
  #              predictorMatrix, formulas, blots, post, c(from, to), 
  #              printFlag)
  ### # INTO THE SAMPLER
  
  ### GO DOWN TO THE SAMPLER FUNCTION ###
  
  ### # AFTER THE SAMPLER
  
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
  
  
## The sampler function:
# sampler <- function(data, m, where, imp, blocks, method, visitSequence, 
#                     predictorMatrix, formulas, blots, post, 
#                     fromto, printFlag, ...)
#   # The sampler controls the actual Gibbs sampling iteration scheme.
#   # This function is called by mice and mice.mids
# {
  fromto <- c(from, to)
  
  from <- fromto[1]
  to <- fromto[2]
  maxit <- to - from + 1
  r <- !is.na(data) # TURE = observed, FALSE = NA
  
  # set up array for convergence checking
  chainMean <- chainVar <- initialize.chain(blocks, maxit, m)
  
  ## THE MAIN LOOP: GIBBS SAMPLER ##
  #######
  # Working parapeters for loops
      k <- 1 # iteration in main loop
      i <- 1 # for each iteration we do things for all the final m datasets
      h <- 1 # for all the variables that should be visited we do things
      j <- blocks[[h]]
  #######    
  
  if (maxit < 1) iteration <- 0
  if (maxit >= 1) {
    if (printFlag)
      cat("\n iter imp variable")
    for (k in from:to) {                   # for the maximum number of iterations
      # begin k loop : main iteration loop
      iteration <- k
      for (i in seq_len(m)) {              # for all of the m datasets that you want to have in the end
        # begin i loop: repeated imputation loop
        if (printFlag)
          cat("\n ", iteration, " ", i)
        
        # prepare the i'th imputation
        # do not overwrite any observed data (fills in missing values)
        for (h in visitSequence) {
          for (j in blocks[[h]]) {
            y <- data[, j]
            ry <- r[, j]
            wy <- where[, j] # locations of missings (True = NA)
            data[(!ry) & wy, j] <- imp[[j]][(!ry)[wy], i]
          }
        }
        
        # impute block-by-block
        for (h in visitSequence) {
          ct <- attr(blocks, "calltype")
          calltype <- ifelse(length(ct) == 1, ct[1], ct[h])
          
          b <- blocks[[h]]
          if (calltype == "formula") ff <- formulas[[h]] else ff <- NULL
          if (calltype == "type") type <- predictorMatrix[h, ] else type <- NULL
          
          user <- blots[[h]]
          
          # univariate/multivariate logic
          theMethod <- method[h]
          empt <- theMethod == ""
          univ <- !empt && !is.passive(theMethod) && 
            !handles.format(paste0("mice.impute.", theMethod))
          mult <- !empt && !is.passive(theMethod) && 
            handles.format(paste0("mice.impute.", theMethod))
          pass <- !empt && is.passive(theMethod) && length(blocks[[h]]) == 1
          if (printFlag & !empt) cat(" ", b)
          
          ## store current state
          oldstate <- get("state", pos = parent.frame())
          newstate <- list(it = k, im = i, 
                           dep = h, 
                           meth = theMethod, 
                           log = oldstate$log)
          assign("state", newstate, pos = parent.frame(), inherits = TRUE)
          
          # (repeated) univariate imputation - type method
          if (univ) {
            for (j in b) {
              imp[[j]][, i] <- 
                sampler.univ(data = data, r = r, where = where, 
                             type = type, formula = ff, 
                             method = theMethod, 
                             yname = j, k = k, 
                             calltype = calltype, 
                             user = user, ...)
              
              data[(!r[, j]) & where[, j], j] <- 
                imp[[j]][(!r[, j])[where[, j]], i]
              
              # optional post-processing
              cmd <- post[j]
              if (cmd != "") {
                eval(parse(text = cmd))
                data[where[, j], j] <- imp[[j]][, i]
              }
            }
          }
          
          # multivariate imputation - type and formula
          if (mult) {
            mis <- !r
            mis[, setdiff(colnames(data), b)] <- FALSE
            data[mis] <- NA
            
            fm <- paste("mice.impute", theMethod, sep = ".")
            if (calltype == "formula")
              imputes <- do.call(fm, args = list(data = data, 
                                                 formula = ff, ...))
            else if (calltype == "type")
              imputes <- do.call(fm, args = list(data = data, 
                                                 type = type, ...))
            else stop("Cannot call function of type ", calltype, 
                      call. = FALSE)
            if (is.null(imputes)) stop("No imputations from ", theMethod, 
                                       h, call. = FALSE)
            for (j in names(imputes)) {
              imp[[j]][, i] <- imputes[[j]]
              data[!r[, j], j] <- imp[[j]][, i]
            }
          }
          
          # passive imputation
          if (pass) {
            for (j in b) {
              wy <- where[, j]
              ry <- r[, j]
              imp[[j]][, i] <- model.frame(as.formula(theMethod), data[wy, ], 
                                           na.action = na.pass)
              data[(!ry) & wy, j] <- imp[[j]][(!ry)[wy], i]
            }
          }
        } # end h loop (blocks)
      }  # end i loop (imputation number)
      
      # store means and sd of m imputes
      k2 <- k - from + 1L
      if (length(visitSequence) > 0L) {
        for (h in visitSequence) {
          for (j in blocks[[h]]) {
            if (!is.factor(data[, j])) {
              chainVar[j, k2, ] <- apply(imp[[j]], 2L, var, na.rm = TRUE)
              chainMean[j, k2, ] <- colMeans(as.matrix(imp[[j]]), na.rm = TRUE)
            }
            if (is.factor(data[, j])) {
              for (mm in seq_len(m)) {
                nc <- as.integer(factor(imp[[j]][, mm], levels = levels(data[, j])))
                chainVar[j, k2, mm] <- var(nc, na.rm = TRUE)
                chainMean[j, k2, mm] <- mean(nc, na.rm = TRUE)
              }
            }
          }
        }
      }
    }  # end main iteration
    
    if (printFlag) {
      r <- get("loggedEvents", parent.frame(1))
      ridge.used <- any(grepl("A ridge penalty", r$out)) 
      if (ridge.used) {
        cat("\n * Please inspect the loggedEvents \n")
      } else {
        cat("\n")
      }
    }
  }
  return(list(iteration = maxit, imp = imp, chainMean = chainMean, chainVar = chainVar))
# }


sampler.univ <- function(data, r, where, type, formula, method, yname, k, 
                         calltype = "type", user, ...) {
  
  # Argumetns
  yname = blocks[[h]]
  # Function
  j <- yname[1L]
  
  if (calltype == "type") {
    vars <- colnames(data)[type != 0]
    pred <- setdiff(vars, j)
    if (length(pred) > 0L) {
      formula <- reformulate(pred, response = j)
      formula <- update(formula, ". ~ . ")
    } else 
      formula <- as.formula(paste0(j, " ~ 1"))
  }
  
  if (calltype == "formula") {
    # move terms other than j from lhs to rhs
    ymove <- setdiff(lhs(formula), j)
    formula <- update(formula, paste(j, " ~ . "))
    if (length(ymove) > 0L) 
      formula <- update(formula, paste("~ . + ", paste(ymove, collapse = "+")))
  }
  
  # get the model matrix
  x <- obtain.design(data, formula) # x is complete (initialized or previous imp)
  
  # expand type vector to model matrix, remove intercept
  if (calltype == "type") {
    type <- type[labels(terms(formula))][attr(x, "assign")]
    x <- x[, -1L, drop = FALSE]
    names(type) <- colnames(x)
  }
  if (calltype == "formula") {
    x <- x[, -1L, drop = FALSE]
    type <- rep(1L, length = ncol(x))
    names(type) <- colnames(x)
  }
  
  # define y, ry and wy
  y <- data[, j] # complete data
  ry <- complete.cases(x, y) & r[, j]
  wy <- complete.cases(x) & where[, j]
  
  # nothing to impute
  if (all(!wy)) return(numeric(0))
  
  cc <- wy[where[, j]]
  if (k == 1L) check.df(x, y, ry)
  
  # remove linear dependencies
  keep <- remove.lindep(x, y, ry)
  x <- x[, keep, drop = FALSE]
  type <- type[keep]
  if (ncol(x) != length(type))
    stop("Internal error: length(type) != number of predictors")
  
  # here we go
  f <- paste("mice.impute", method, sep = ".")
  imputes <- data[wy, j]
  imputes[!cc] <- NA
  
  args <- c(list(y = y, ry = ry, x = x, wy = wy, type = type), user)
  imputes[cc] <- do.call(f, args = args) # f is the name of the imputation function to be called. 
  # hence duplicate the file for mive.impute.norm.boot, and substitute the algorithm from the CALIBRE package
  # to obtain the same algorithm with a new method
  imputes
}
