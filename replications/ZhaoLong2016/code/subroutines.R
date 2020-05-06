### Title:    Replication Zhao Long 2016 - Functions
### Author:   Edoardo Costantini
### Created:  2020-02-20
### Modified: 2020-02-20

# Functions ---------------------------------------------------------------

## Run one replication of the simulation:
doRep <- function(rp, conds, parms) { 
  ## Setup the PRNG: (for setting seed with guaranteed independece)
  # form the rlecuyer package
  # rp is the replication index, in the mcapply function it the argument
  # that takes the elements of the vector specified for X
  # This creates names for the streams so that you can recover them
  .lec.SetPackageSeed(rep(parms$seed, 6))
  .lec.CreateStream(c(1 : parms$nStreams)) # create 1000 streams
  .lec.CurrentStream(rp) # use the rp sequence out of the 1000
  
  ## Progress report - Start
  parms$rep_counter <- parms$rep_counter + 1 # increase progres report counter
  cat(paste0(Sys.time(), " - Starts Repetition: ", rp, 
             "\n"),
      file = paste0(parms$outDir, parms$report_file_name),
      append = TRUE)
  
  ## Do the calculations for each set of crossed conditions:
  rp_out <- vector("list", nrow(conds)) # store output of repetition
    names(rp_out) <- paste0("cond_", paste0(conds[, 1], "_",  conds[, 3]) )
    
  for(i in 1 : nrow(conds)) {
    # Perform runCell until you get no error
    # there are some cases where the dataset that you have generated
    # has a IURR lasso that selects more variables than cases. In that
    # situation you need to discard the cycle and repeat it. The while
    # does just that! Not elegant but should work for now.
    
    exit_while <- "no"
    
    while(exit_while == "no"){
      
      rp_out[[i]] <- try(runCell(cond = conds[i, ], 
                              m = parms$chains, 
                              iters = parms$iters), 
                      silent = TRUE)
      
      if (class(rp_out[[i]]) != "try-error") {
        exit_while <- "yes"
      }
      
    }
    
  }
  ## Return Function output  
  return(rp_out)
}

runCell <- function(cond, m = 5, iters = 1) {
  ## Description
  # Generates data for 1 condition and performs imutations according to
  # DURR, IURR, BLasso, and basic methods.
  ## For internals
  # source("./init.R")
  # source("./functions.R")
  # source("./fun_DURR_impute.R")
  # source("./fun_IURR_impute.R")
  # source("./fun_BLasso_impute.R")
  # source("./fun_blassoHans_impute.R")
  # m = parms$chains
  # iters = 1
  # cond <- conds[1, ]
  
  ## Data ------------------------------------------------------------------ ##
  # Gen 1 dataset w/ missing values; 1 fully-obs data for out-of-sample rmse
  
  Xy <- genData(cond, parms)
  Xy_mis <- imposeMiss(Xy, parms)$Xy_miss
  
  miss_descrps <- mean(is.na(Xy_mis[, "z1"]))
  
  ## Imputation ------------------------------------------------------------ ##
  # Impute m times the data w/ missing values w/ different methods
  
  # Impute according to DURR method
  imp_DURR_lasso <- impute_DURR(Xy_mis = Xy_mis, 
                                cond = cond, 
                                chains = m, 
                                iters = iters, 
                                reg_type="lasso")
  
  # Impute according to IURR method
  imp_IURR_lasso <- impute_IURR(Xy_mis = Xy_mis, 
                                cond = cond, 
                                chains = m, 
                                reg_type = "lasso")
  
  # Impute according to Hans Blasso method
  imp_BLasso <- impute_BLAS_hans(Xy = Xy, Xy_mis = Xy_mis, 
                                 chains = m, 
                                 iter_bl = parms$iter_bl, 
                                 burn_bl = parms$burn_bl)
  
  # MICE w/ true model
  S <- parms$S_all[[ which(paste0("q", cond[3]) == names(parms$S_all)) ]]
  varTRUE <- c(1, (S+1), ncol(Xy_mis))
  imp_MI_T_mids <- mice::mice(Xy_mis[, varTRUE], 
                                # not elegant way of selecting y and active set
                              m = m,
                              maxit = iters,
                              method = "norm")
  imp_MI_T <- mice::complete(imp_MI_T_mids, "all")
  
  # MICE-50 w/ ridge prior
  var50 <- get_50_best(Xy_mis = Xy_mis, S = S) # gets50 most corr w/ z1 (+ z1)
  imp_MI_50_mids <- mice::mice(Xy_mis[, var50],
                               m = m,
                               maxit = iters,
                               ridge = 1e-5,
                               method = "norm")
  imp_MI_50 <- mice::complete(imp_MI_50_mids, "all")
  
  ## Analyse --------------------------------------------------------------- ##
  # For each imp method, analyse all datasets based on model defined in init.R
  fits_md <- lapply(list(imp_DURR_lasso,
                         imp_IURR_lasso,
                         imp_BLasso,
                         imp_MI_T,
                         imp_MI_50), 
                    fit_models, mod = parms$formula)
  
  # CC (complete case analysis)
  lm_fit_cc <- lm(parms$formula, data = Xy_mis, na.action = na.omit)
  bs_cc <- coef(lm_fit_cc)
  
  ## Pooling --------------------------------------------------------------- ##
  # For each imp method, pool estimates across the m datasets
  pool_EST <- lapply(fits_md, get_pool_EST)
  pool_CI  <- lapply(fits_md, get_pool_CI)
  
  # append CC results
  pool_EST[[length(pool_EST)+1]] <- coef(lm_fit_cc)
  pool_CI[[length(pool_CI)+1]] <- confint(lm_fit_cc)
  
  # give meaningful names
  names(pool_EST) <- parms$methods
  names(pool_CI) <- parms$methods
  
  # compute bias and check coverage
  cond_bias <- sapply(pool_EST, bias_est, x_true = parms$b)
  cond_CIco <- sapply(pool_CI, check_cover)
  
  ## Store output ---------------------------------------------------------- ##
  output <- list(pool_est = pool_EST,
                 pool_conf = pool_CI,
                 miss_descrps = miss_descrps,
                 cond_bias = cond_bias,
                 cond_CIco = cond_CIco,
                 parms = parms)
  
  return(output)
}
