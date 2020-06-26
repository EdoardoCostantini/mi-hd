### Title:    Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-05-19

# Functions ---------------------------------------------------------------

## Run one replication of the simulation:
doRep <- function(rp, conds, parms) {
  ## For 1 repetition, performs simulation across all conditions
  ## Setup the PRNG: (for setting seed with guaranteed independece)
  # form the rlecuyer package
  # rp is the replication index, in the mcapply function it the argument
  # that takes the elements of the vector specified for X
  # This creates names for the streams so that you can recover them
  
  # .lec.SetPackageSeed(rep(parms$seed, 6))
  # .lec.CreateStream(c(1 : parms$nStreams)) # create 1000 streams
  # .lec.CurrentStream(rp) # use the rp sequence out of the 1000
 
  ## Progress report - Start
  # parms$rep_counter <- parms$rep_counter + 1 # increase progres report counter
  cat(paste0(Sys.time(), " - Starts Repetition: ", rp, 
             "\n"),
      file = paste0(parms$outDir, parms$report_file_name),
      append = TRUE)
  
  ## Do the calculations for each set of crossed conditions:
  rp_out <- vector("list", nrow(conds)) # store output of repetition
    names(rp_out) <- paste0("cond_", paste0(conds[, "p"], "_",  conds[, "pm"]) )
    
  for(i in 1 : nrow(conds)) {
    # Perform runCell until you get no error
    # there are some cases where the dataset that you have generated
    # has a IURR lasso that selects more variables than cases. In that
    # situation you need to discard the cycle and repeat it. The while
    # does just that! Not elegant but should work for now.
    
    exit_while <- "no"
    
    while(exit_while == "no"){

      rp_out[[i]] <- tryCatch(
        {
          # Try running simulation for condition i, repetition rp
          runCell(cond = conds[i, ], 
                  parms = parms,
                  rep_status = rp)
        }, 
        error = function(report) {
          err <- paste0("Original Error: ", report)
          return(err)
        }
      )
      
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
  # Given 1 condition, Generates 1 dataset and performs imutations according to
  # selected methods
  ## For internals
  # source("./init.R")
  # cond <- conds[1, ]; # set.seed(2020)
  
  ## Data ------------------------------------------------------------------ ##
  # Gen 1 dataset w/ missing values; 1 fully-obs data for out-of-sample rmse
  if(cond$latent == FALSE){
    Xy <- genDt_mvn(cond, parms)
    Xy_mis <- imposeMiss(Xy, parms, cond)
    Xy_mis <- cbind(Xy_mis[parms$z_m_id], Xy_mis[-parms$z_m_id])
    
    O <- !is.na(Xy_mis) # matrix index of observed values
    miss_descrps <- colMeans(!O[, 1:parms$zm_n]) 
      # check correct proportion of missing cases
  } else {
    Xy <- genDataCat(cond, parms)
    Xy_mis <- imposeMissCat(Xy, parms)$Xy_miss
    
    O <- !is.na(Xy_mis) # matrix index of observed values
    miss_descrps <- mean(rowSums(!O) != 0) # check correct miss %
  }
  
  ## Imputation ------------------------------------------------------------ ##
  # Impute m times the data w/ missing values w/ different methods

  # Impute according to DURR method
  # Ridge
  imp_DURR_rd <- impute_DURR(Z = Xy_mis,
                             O = as.data.frame(O),
                             reg_type = "ridge",
                             perform = parms$meth_sel$DURR_rd,
                             parms = parms)
  update_report("DURR ridge", rep_status, parms)
  
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
                                 parms = parms,
                                 perform = parms$meth_sel$blasso)
  update_report("blasso", rep_status, parms)

  # Impute according to Howard Et Al 2015 PCA appraoch
  sem_cnvg <- FALSE
  while (sem_cnvg == FALSE) {
    imp_PCA <- impute_PCA(Z = Xy_mis, parms = parms)
    imp_PCA_fit <- tryCatch({fit_sat_model(imp_PCA$dats)}, 
                            error = function(report) {
                              err <- paste0("Original Error: ", report)
                              return(err)
                            },
                            warning = function(report) {
                              err <- paste0("Original Warning: ", report)
                              return(err)
                            })
    sem_cnvg <- !is.character(imp_PCA_fit)
  }
  update_report("MICE-PCA", rep_status, parms)
  
  # MICE-CART traditional
  imp_CART <- impute_CART(Z = Xy_mis,
                          O = as.data.frame(O),
                          perform = parms$meth_sel$MI_CART,
                          parms = parms)
  
  update_report("MICE-CART", rep_status, parms)
  
  # MICE-RF
  imp_RANF <- impute_RANF(Z = Xy_mis,
                          O = as.data.frame(O),
                          perform = parms$meth_sel$MI_RF,
                          parms = parms)
  
  update_report("Random Forest", rep_status, parms)
  
  # MICE w/ true model
  sem_cnvg <- FALSE
  while (sem_cnvg == FALSE) {
    imp_MICE_TR <- impute_MICE_TR(Z = Xy_mis,
                                  cond = cond,
                                  parms = parms)
    imp_MICE_TR_fit <- tryCatch({fit_sat_model(imp_MICE_TR$dats)}, 
                            error = function(report) {
                              err <- paste0("Original Error: ", report)
                              return(err)
                            },
                            warning = function(report) {
                              err <- paste0("Original Warning: ", report)
                              return(err)
                            })
    sem_cnvg <- !is.character(imp_PCA_fit)
  }
  
  update_report("MICE-TR", rep_status, parms)
  
  # missForest
  imp_missFor <- impute_missFor(Z = Xy_mis, parms = parms)
  
  ## Convergence ----------------------------------------------------------- ##
  
  imp_values <- list(DURR_rd = imp_DURR_rd$imps,
                     DURR_la = imp_DURR_la$imps,
                     DURR_el = imp_DURR_el$imps,
                     IURR_la = imp_IURR_la$imps,
                     IURR_el = imp_IURR_el$imps,
                     blasso = imp_blasso$imps,
                     MICE_PCA = imp_PCA$mids,
                     MICE_CART = imp_CART$mids,
                     MICE_RF = imp_RANF$imps,
                     MICE_TR = imp_MICE_TR$mids) # I need this for convergence
  
  ## Analyse --------------------------------------------------------------- ##
  # For each imp method, analyse all datasets based on model defined in init.R
  fits_md <- lapply(list(imp_DURR_rd$dats,
                         imp_DURR_la$dats,
                         imp_DURR_el$dats,
                         imp_IURR_la$dats,
                         imp_IURR_el$dats,
                         imp_blasso$dats,
                         imp_PCA$dats,
                         imp_CART$dats,
                         imp_RANF$dats,
                         imp_MICE_TR$dats), 
                    fit_sat_model)
  
  # missForest single imputation
  fit_mf <- sem(parms$lav_model, 
                data = imp_missFor$dats, 
                likelihood = "wishart")
  
  # Gold Standard
  fit_gs <- sem(parms$lav_model, 
                data = Xy, 
                likelihood = "wishart")
  
  # CC (complete case analysis)
  fit_cc <- sem(parms$lav_model, 
                data = Xy_mis[rowSums(!O) == 0, ], # select complete casese
                likelihood = "wishart")
  
  ## Pooling --------------------------------------------------------------- ##
  # For each imp method, pool estimates across the m datasets
  pool_EST <- sapply(fits_md[lapply(fits_md, length) != 0], 
                     get_pool_EST)
  pool_CI  <- sapply(fits_md[lapply(fits_md, length) != 0], 
                     get_pool_CI)
  
  # append single imputations, and GS and CC results
  all_EST <- cbind(pool_EST, 
                get_EST(list(fit_mf, fit_gs, fit_cc)))
  
  all_CI <- cbind(pool_CI, 
               get_CI(list(fit_mf, fit_gs, fit_cc)))
  
  # give meaningful names
  colnames(all_EST) <- parms$methods
  colnames(all_CI) <- parms$methods

  # aggregate times
  imp_time <- list(DURR_rd   = imp_DURR_rd$time,
                   DURR_la   = imp_DURR_la$time,
                   DURR_el   = imp_DURR_el$time,
                   IURR      = imp_IURR_la$time,
                   IURR_el   = imp_IURR_el$time,
                   blasso    = imp_blasso$time,
                   MICE_PCA  = imp_PCA$time,
                   MICE_CART = imp_CART$time,
                   MICE_RF   = imp_RANF$time,
                   MICE_TR   = imp_MICE_TR$time)
  imp_time <- do.call(cbind, imp_time)[1,]
  
  ## Store output ---------------------------------------------------------- ##
  output <- list(cond = cond,
                 dat_full = Xy,
                 dat_miss = Xy_mis,
                 all_EST = all_EST,
                 all_CI = all_CI,
                 miss_descrps = miss_descrps,
                 run_time_min = imp_time,
                 imp_values = imp_values)
  return(output)
}
