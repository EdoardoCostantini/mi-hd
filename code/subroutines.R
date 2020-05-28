### Title:    Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-05-19

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
  # parms$rep_counter <- parms$rep_counter + 1 # increase progres report counter
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
                                 parms = parms,
                                 rep_status = rp), 
                      silent = TRUE)
      
      if (class(rp_out[[i]]) != "try-error") {
        exit_while <- "yes"
      }
      
    }
    
  }
  ## Return Function output  
  return(rp_out)
}

runCell <- function(cond, parms, rep_status) {
  ## Description
  # Generates data for 1 condition and performs imutations according to
  # DURR, IURR, BLasso, and basic methods.
  ## For internals
  # source("./init.R")
  # cond <- conds[1, ]; set.seed(1234)
  
  ## Data ------------------------------------------------------------------ ##
  # Gen 1 dataset w/ missing values; 1 fully-obs data for out-of-sample rmse
  Xy <- genData(cond, parms)
  Xy_mis <- imposeMiss(Xy, parms)$Xy_miss
  
  O <- !is.na(Xy_mis) # matrix index of observed values
  miss_descrps <- mean(rowSums(!O) != 0) # check correct miss %
  
  ## Imputation ------------------------------------------------------------ ##
  # Impute m times the data w/ missing values w/ different methods

  # Impute according to DURR method
    # Lasso
  imp_DURR_la <- impute_DURR(Z = Xy_mis,
                             O = as.data.frame(O),
                             reg_type = "lasso",
                             perform = parms$meth_sel$DURR_la,
                             parms = parms)
  update_report("DURR lasso", rep_status, parms)
  
    # Elastic Net
  imp_DURR_el <- impute_DURR(Z = Xy_mis,
                             O = as.data.frame(O),
                             reg_type = "el",
                             perform = parms$meth_sel$DURR_el,
                             parms = parms)
  update_report("DURR elastic net", rep_status, parms)

  # Impute according to IURR method
    # Lasso
  imp_IURR_la <- impute_IURR(Z = Xy_mis,
                                O = as.data.frame(O),
                                reg_type = "lasso",
                                perform = parms$meth_sel$IURR_la,
                                parms = parms)
  update_report("IURR lasso", rep_status, parms)
  
    # Elastic Net
  imp_IURR_el <- impute_IURR(Z = Xy_mis,
                             O = as.data.frame(O),
                             reg_type = "el",
                             perform = parms$meth_sel$IURR_el,
                             parms = parms)
  update_report("IURR elastic net", rep_status, parms)

  # Impute according to Hans Blasso method
  imp_blasso <- impute_BLAS_hans(Z = Xy_mis, 
                                 O = as.data.frame(O),
                                 parms = parms)
  update_report("blasso", rep_status, parms)

  # Impute according to Howard Et Al 2015 PCA appraoch
  imp_PCA <- impute_PCA(Z = Xy_mis, parms = parms)
  update_report("MICE-PCA", rep_status, parms)
  
  # MICE-CART traditional
  imp_MICE_CART <- impute_MICE_CART(Z = Xy_mis,
                                    parms = parms)
  
  update_report("MICE-CART", rep_status, parms)
  
  # MICE-CART bootstarp sample
  imp_MICE_CARTbb <- impute_MICE_CARTbb(Z = Xy_mis,
                                        parms = parms)

  update_report("MICE-CARTbb", rep_status, parms)
  
  # MICE-RF
  imp_MICE_RF <- impute_MICE_RF(Z = Xy_mis,
                                parms = parms)

  update_report("MICE-RF", rep_status, parms)
  
  # MICE w/ true model
  imp_MICE_TR <- impute_MICE_TR(Z = Xy_mis,
                                AS_size = cond[3],
                                parms = parms)
  
  update_report("MICE-TR", rep_status, parms)
  
  ## Convergence ----------------------------------------------------------- ##
  
  imp_values <- list(DURR = imp_DURR_la$imps,
                     DURR_el = imp_DURR_el$imps,
                     IURR = imp_IURR_la$imps,
                     IURR_el = imp_IURR_el$imps,
                     blasso = imp_blasso$imps,
                     MICE_PCA = imp_PCA$imps,
                     MICE_CART = imp_MICE_CART$imps,
                     MICE_CARTbb = imp_MICE_CARTbb$imps,
                     MICE_RF = imp_MICE_RF$imps,
                     MICE_TR = imp_MICE_TR$mids)
  
  ## Analyse --------------------------------------------------------------- ##
  # For each imp method, analyse all datasets based on model defined in init.R
  fits_md <- lapply(list(imp_DURR_la$dats,
                         imp_DURR_el$dats,
                         imp_IURR_la$dats,
                         imp_IURR_el$dats,
                         imp_blasso$dats,
                         imp_PCA$dats,
                         imp_MICE_CART$dats,
                         imp_MICE_CARTbb$dats,
                         imp_MICE_RF$dats,
                         imp_MICE_TR$dats), 
                    fit_models, mod = parms$formula)
  
  # Gold Standard
  lm_fit_gs <- lm(parms$formula, data = Xy)
  bs_gs <- coef(lm_fit_gs)
  CI_gs <- confint(lm_fit_gs)
  
  # CC (complete case analysis)
  lm_fit_cc <- lm(parms$formula, data = Xy_mis, na.action = na.omit)
  bs_cc <- coef(lm_fit_cc)
  CI_cc <- confint(lm_fit_cc)
  
  ## Pooling --------------------------------------------------------------- ##
  # For each imp method, pool estimates across the m datasets
  pool_EST <- lapply(fits_md[lapply(fits_md, length) != 0], 
                     get_pool_EST)
  pool_CI  <- lapply(fits_md[lapply(fits_md, length) != 0], 
                     get_pool_CI)
  
  # append GS and CC results
  pool_EST <- c(pool_EST, list(bs_gs, bs_cc))
  pool_CI <- c(pool_CI, list(CI_gs, CI_cc))
  
  # give meaningful names
  names(pool_EST) <- parms$methods
  names(pool_CI) <- parms$methods
  
  # compute bias and check coverage
  cond_bias <- sapply(pool_EST, bias_est, x_true = parms$b)
  cond_CIco <- sapply(pool_CI, check_cover)
  
  # aggregate times
  imp_time <- list(DURR        = imp_DURR_la$time,
                   DURR_el     = imp_DURR_el$time,
                   IURR        = imp_IURR_la$time,
                   IURR_el     = imp_IURR_el$time,
                   blasso      = imp_blasso$time,
                   MICE_PCA    = imp_PCA$time,
                   MICE_CART   = imp_MICE_CART$time,
                   MICE_CARTbb = imp_MICE_CARTbb$time,
                   MICE_RF     = imp_MICE_RF$time,
                   MICE_TR     = imp_MICE_TR$time)
  imp_time <- do.call(cbind, imp_time)[1,]
  
  ## Store output ---------------------------------------------------------- ##
  output <- list(pool_est = pool_EST,
                 pool_conf = pool_CI,
                 miss_descrps = miss_descrps,
                 cond_bias = cond_bias,
                 cond_CIco = cond_CIco,
                 run_time_min = imp_time,
                 imp_values = imp_values, 
                 cond = cond,
                 parms = parms)
  return(output)
}
