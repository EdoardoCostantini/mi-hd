### Title:    Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-05-19

# Functions ---------------------------------------------------------------

## Run one replication of the simulation:
doRep <- function(rp, conds, parms, debug = FALSE) {
  ## For internals
  # rp = 1
  # debug = FALSE
  
  ## For 1 repetition, performs simulation across all conditions
  ## Setup the PRNG: (for setting seed with guaranteed independece)
  # form the rlecuyer package
  # rp is the replication index, in the mcapply function it the argument
  # that takes the elements of the vector specified for X
  # This creates names for the streams so that you can recover them
  
  .lec.SetPackageSeed(rep(parms$seed, 6))
  if(!rp %in% .lec.GetStreams()) # if the streams do not exist yet
    .lec.CreateStream(c(1 : parms$nStreams)) # then
  .lec.CurrentStream(rp) # this is equivalent to setting the seed Rle
 
  ## Progress report - Start
  # parms$rep_counter <- parms$rep_counter + 1 # increase progres report counter
  cat(paste0(Sys.time(), " - Starts Repetition: ", rp, 
             "\n"),
      file = paste0(parms$outDir, parms$report_file_name),
      append = TRUE)
  
  ## Do the calculations for each set of crossed conditions:
  rp_out <- vector("list", nrow(conds)) # store output of repetition
    # names(rp_out) <- paste0("cond_", paste0(conds[, "p"], "_",  conds[, "pm"]) )
    names(rp_out) <- paste0("cond_", 
                            sapply(1:nrow(conds), 
                                   function(x) paste0(conds[x,], 
                                                      collapse = "_")
                                   )
                            )
    
  for(i in 1 : nrow(conds)) {
    
  if(debug == TRUE){
    rp_out[[i]] <- capture.output(tryCatch({
      withErrorTracing({runCell(cond = conds[i, ],
                                parms = parms,
                                rep_status = rp)})
    }, error = function(e){
      e <<- e
      cat("ERROR: ", e$message, "\nin ")
      print(e$call)
    }))
  } else {
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
  }
    
  }
  ## Return Function output  
  return(rp_out)
}

# Experiment 1 ------------------------------------------------------------

runCell <- function(cond, parms, rep_status) {
  ## Description
  # Given 1 condition, Generates 1 dataset and performs imutations according to
  # selected methods
  ## For internals
  # source("./init.R")
  # cond <- conds[1, ]
  
  ## Data ------------------------------------------------------------------ ##
  # Gen 1 dataset w/ missing values; 1 fully-obs data for out-of-sample rmse
  Xy <- simData_exp1(cond, parms)
  Xy_mis <- imposeMiss(Xy, parms, cond)
  Xy_mis <- cbind(Xy_mis[, parms$z_m_id], Xy_mis[, -parms$z_m_id])
  
  O <- !is.na(Xy_mis) # matrix index of observed values
  miss_descrps <- colMeans(!O[, 1:parms$zm_n]) 
  
  ## Imputation ------------------------------------------------------------ ##
  # Impute m times the data w/ missing values w/ different methods

  # Impute according to DURR method
  
    # Lasso
  imp_DURR_la <- impute_DURR(Z = Xy_mis,
                             O = as.data.frame(O),
                             reg_type = "lasso",
                             perform = parms$meth_sel$DURR_la,
                             parms = parms)

  update_report("DURR lasso", rep_status, parms, 
                cnd = cond,
                perform = parms$meth_sel$DURR_la)
  
    # Elastic Net
  imp_DURR_el <- impute_DURR(Z = Xy_mis,
                             O = as.data.frame(O),
                             reg_type = "el",
                             perform = parms$meth_sel$DURR_el,
                             parms = parms)
  update_report("DURR elastic net", rep_status, parms, 
                cnd = cond,
                perform = parms$meth_sel$DURR_el)

  # Impute according to IURR method
    # Lasso
  imp_IURR_la <- impute_IURR(Z = Xy_mis,
                                O = as.data.frame(O),
                                reg_type = "lasso",
                                perform = parms$meth_sel$IURR_la,
                                parms = parms)
  update_report("IURR lasso", rep_status, parms, 
                cnd = cond,
                perform = parms$meth_sel$IURR_la)
  
    # Elastic Net
  imp_IURR_el <- impute_IURR(Z = Xy_mis,
                             O = as.data.frame(O),
                             reg_type = "el",
                             perform = parms$meth_sel$IURR_el,
                             parms = parms)
  update_report("IURR elastic net", rep_status, parms, 
                cnd = cond,
                perform = parms$meth_sel$IURR_el)

  # Impute according to Hans Blasso method
  imp_blasso <- impute_BLAS_hans(Z = Xy_mis, 
                                 O = as.data.frame(O),
                                 parms = parms,
                                 perform = parms$meth_sel$blasso)
  update_report("blasso", rep_status, parms, 
                cnd = cond,
                perform = parms$meth_sel$blasso)
  
  # Impute according to van Buuren Ridge
  imp_bridge <- impute_BRIDGE(Z = Xy_mis, 
                              O = as.data.frame(O),
                              parms = parms,
                              perform = parms$meth_sel$bridge)
  update_report("bridge", rep_status, parms, 
                cnd = cond,
                perform = parms$meth_sel$bridge)

  # Impute according to Howard Et Al 2015 PCA appraoch
  imp_PCA <- impute_PCA(Z = Xy_mis, parms = parms)
  update_report("MICE-PCA", rep_status, parms, 
                cnd = cond,
                perform = parms$meth_sel$MI_PCA)
  
  # MICE-CART traditional
  imp_CART <- impute_CART(Z = Xy_mis,
                          O = as.data.frame(O),
                          perform = parms$meth_sel$MI_CART,
                          parms = parms)
  
  update_report("MICE-CART", rep_status, parms, 
                cnd = cond,
                perform = parms$meth_sel$MI_CART)
  
  # MICE-RF
  imp_RANF <- impute_RANF(Z = Xy_mis,
                          O = as.data.frame(O),
                          perform = parms$meth_sel$MI_RF,
                          parms = parms)
  
  update_report("Random Forest", rep_status, parms, 
                cnd = cond,
                perform = parms$meth_sel$MI_RF)
  
  # MICE w/ true model
  imp_MICE_OP <- impute_MICE_OP(Z = Xy_mis,
                                cond = cond,
                                perform = parms$meth_sel$MI_OP,
                                parms = parms)
  update_report("MICE-TR", rep_status, parms, 
                cnd = cond,
                perform = parms$meth_sel$MI_OP)
  
  # missForest
  imp_missFor <- impute_missFor(Z = Xy_mis, parms = parms)
  missForest_out <- capture.output(impute_missFor(Z = Xy_mis, parms = parms))
  missFprest_outdt <- capture.output(imp_missFor$dats)
  
  ## Convergence ----------------------------------------------------------- ##
  
  imp_values <- list(DURR_la = imp_DURR_la$imps,
                     DURR_el = imp_DURR_el$imps,
                     IURR_la = imp_IURR_la$imps,
                     IURR_el = imp_IURR_el$imps,
                     bridge  = imp_bridge$imps,
                     blasso  = imp_blasso$imps,
                     MI_PCA  = imp_PCA$mids,
                     MI_CART = imp_CART$imps,
                     MI_RF   = imp_RANF$imps,
                     MI_OP   = imp_MICE_OP$mids) # I need this for convergence
  
  ## Analyse --------------------------------------------------------------- ##
  # For each imp method, analyse all datasets based on model defined in init.R
  
  ## MLE estaimte of mean, variance, covariance
  
  # Multiple imputed dataset
  
  sem_fits <- lapply(list(DURR_la = imp_DURR_la$dats,
                          DURR_el = imp_DURR_el$dats,
                          IURR_la = imp_IURR_la$dats,
                          IURR_el = imp_IURR_el$dats,
                          bridge  = imp_bridge$dats,
                          blasso  = imp_blasso$dats,
                          MI_PCA  = imp_PCA$dats,
                          MI_CART = imp_CART$dats,
                          MI_RF   = imp_RANF$dats,
                          MI_OP   = imp_MICE_OP$dats), 
                     fit_sat_model)
  
  # Single dataset
  sem_sndt <- lapply(list(missFor = imp_missFor$dats,          
                          GS      = Xy,                        
                          CC      = Xy_mis[rowSums(!O) == 0, ]), 
                     sem, model = parms$lav_model, likelihood = "wishart")
  
  sem_sndt_output <- capture.output(sem_sndt)
  
  ## LM model
  
  # Multiple datasets
  
  lm_fits <- lapply(list(DURR_la = imp_DURR_la$dats,
                         DURR_el = imp_DURR_el$dats,
                         IURR    = imp_IURR_la$dats,
                         IURR_el = imp_IURR_el$dats,
                         bridge  = imp_bridge$dats,
                         blasso  = imp_blasso$dats,
                         MI_PCA  = imp_PCA$dats,
                         MI_CART = imp_CART$dats,
                         MI_RF   = imp_RANF$dats,
                         MI_OP   = imp_MICE_OP$dats), 
                    fit_lm_models, vrbs = parms$lm_model)
  
  # Single dataset
  
  lm_sndt <- fit_lm_models(list(missFor = imp_missFor$dats, 
                                GS      = Xy, 
                                CC      = Xy_mis[rowSums(!O) == 0, ]),
                           vrbs = parms$lm_model)
  cov
  ## Pooling --------------------------------------------------------------- ##
  # For each imp method, pool estimates across the m datasets
  
  # MLE mean, var, cov
  
  sem_pool_EST <- sapply(sem_fits[lapply(sem_fits, length) != 0], 
                         sem_pool_EST_f)
  sem_pool_CI  <- sapply(sem_fits[lapply(sem_fits, length) != 0], 
                         sem_pool_CI_f)
  
  # append single imputations, and GS and CC results
  sem_pool_EST <- cbind(sem_pool_EST,
                        sem_EST(sem_sndt))
  
  sem_pool_CI <- cbind(sem_pool_CI,
                       sem_CI(sem_sndt))
  
  # LM models
  lm_pool_est <- sapply(lm_fits[lapply(lm_fits, length) != 0], lm_pool_EST_f)
  lm_pool_CI <- sapply(lm_fits[lapply(lm_fits, length) != 0], lm_pool_CI_f)
  
  # append single imputations, and GS and CC results
  lm_pool_EST <- cbind(lm_pool_est,
                       lm_EST(lm_sndt))
  
  lm_pool_CI <- cbind(lm_pool_CI,
                      lm_CI(lm_sndt))
  
  
  ## Times ----------------------------------------------------------------- ##
  # aggregate imputation times
  
  imp_time <- list(DURR_la = imp_DURR_la$time,
                   DURR_el = imp_DURR_el$time,
                   IURR    = imp_IURR_la$time,
                   IURR_el = imp_IURR_el$time,
                   bridge  = imp_bridge$time,
                   blasso  = imp_blasso$time,
                   MI_PCA  = imp_PCA$time,
                   MI_CART = imp_CART$time,
                   MI_RF   = imp_RANF$time,
                   MI_OP   = imp_MICE_OP$time)
  imp_time <- do.call(cbind, imp_time)[1,]
  
  ## Store output ---------------------------------------------------------- ##
  output <- list(cond = cond,
                 dat_full     = Xy,
                 dat_miss     = Xy_mis,
                 sem_EST      = sem_pool_EST,
                 sem_CI       = sem_pool_CI,
                 lm_EST       = lm_pool_EST,
                 lm_CI        = lm_pool_CI,
                 miss_descrps = miss_descrps,
                 run_time_min = imp_time,
                 imp_values   = imp_values)
  return(output)
}


# Experiment 3 ------------------------------------------------------------

runCell_exp3 <- function(cond, parms, rep_status) {
  ## Description
  # Given 1 condition, Generates 1 dataset and performs imutations according to
  # selected methods
  ## For internals
  # source("./init.R")
  # cond <- conds[1, ]
  
  ## Data ------------------------------------------------------------------ ##
  # Gen 1 dataset w/ missing values; 1 fully-obs data for out-of-sample rmse
  simData_list <- simData_exp3(parms, cond)
    Xy <- simData_list$dat
  Xy_mis <- imposeMiss_lv(simData_list, parms, cond)
  
  O <- !is.na(Xy_mis) # matrix index of observed values
  miss_descrps <- colMeans(!O[, 1:parms$zm_n])
  
  ## Imputation ------------------------------------------------------------ ##
  # Impute m times the data w/ missing values w/ different methods
  
  # Impute according to DURR method
  
  # Lasso
  imp_DURR_la <- impute_DURR(Z = Xy_mis,
                             O = as.data.frame(O),
                             reg_type = "lasso",
                             perform = parms$meth_sel$DURR_la,
                             parms = parms)
  
  update_report("DURR lasso", rep_status, parms, 
                cnd = cond,
                perform = parms$meth_sel$DURR_la)
  
  # Elastic Net
  imp_DURR_el <- impute_DURR(Z = Xy_mis,
                             O = as.data.frame(O),
                             reg_type = "el",
                             perform = parms$meth_sel$DURR_el,
                             parms = parms)
  update_report("DURR elastic net", rep_status, parms, 
                cnd = cond,
                perform = parms$meth_sel$DURR_el)
  
  # Impute according to IURR method
  # Lasso
  imp_IURR_la <- impute_IURR(Z = Xy_mis,
                             O = as.data.frame(O),
                             reg_type = "lasso",
                             perform = parms$meth_sel$IURR_la,
                             parms = parms)
  update_report("IURR lasso", rep_status, parms, 
                cnd = cond,
                perform = parms$meth_sel$IURR_la)
  
  # Elastic Net
  imp_IURR_el <- impute_IURR(Z = Xy_mis,
                             O = as.data.frame(O),
                             reg_type = "el",
                             perform = parms$meth_sel$IURR_el,
                             parms = parms)
  update_report("IURR elastic net", rep_status, parms, 
                cnd = cond,
                perform = parms$meth_sel$IURR_el)
  
  # Impute according to Hans Blasso method
  imp_blasso <- impute_BLAS_hans(Z = Xy_mis, 
                                 O = as.data.frame(O),
                                 parms = parms,
                                 perform = parms$meth_sel$blasso)
  update_report("blasso", rep_status, parms, 
                cnd = cond,
                perform = parms$meth_sel$blasso)
  
  # Impute according to van Buuren Ridge
  imp_bridge <- impute_BRIDGE(Z = Xy_mis, 
                              O = as.data.frame(O),
                              parms = parms,
                              perform = parms$meth_sel$bridge)
  update_report("bridge", rep_status, parms, 
                cnd = cond,
                perform = parms$meth_sel$bridge)
  
  # Impute according to Howard Et Al 2015 PCA appraoch
  imp_PCA <- impute_PCA(Z = Xy_mis, parms = parms)
  update_report("MICE-PCA", rep_status, parms, 
                cnd = cond,
                perform = parms$meth_sel$MI_PCA)
  
  # MICE-CART traditional
  imp_CART <- impute_CART(Z = Xy_mis,
                          O = as.data.frame(O),
                          perform = parms$meth_sel$MI_CART,
                          parms = parms)
  
  update_report("MICE-CART", rep_status, parms, 
                cnd = cond,
                perform = parms$meth_sel$MI_CART)
  
  # MICE-RF
  imp_RANF <- impute_RANF(Z = Xy_mis,
                          O = as.data.frame(O),
                          perform = parms$meth_sel$MI_RF,
                          parms = parms)
  
  update_report("Random Forest", rep_status, parms, 
                cnd = cond,
                perform = parms$meth_sel$MI_RF)
  
  # MICE w/ true model
  imp_MICE_OP <- impute_MICE_OP(Z = Xy_mis,
                                cond = cond,
                                perform = parms$meth_sel$MI_OP,
                                parms = parms)
  update_report("MICE-TR", rep_status, parms, 
                cnd = cond,
                perform = parms$meth_sel$MI_OP)
  
  # missForest
  imp_missFor <- impute_missFor(Z = Xy_mis, parms = parms)
  missForest_out <- capture.output(impute_missFor(Z = Xy_mis, parms = parms))
  missFprest_outdt <- capture.output(imp_missFor$dats)
  
  ## Convergence ----------------------------------------------------------- ##
  
  imp_values <- list(DURR_la = imp_DURR_la$imps,
                     DURR_el = imp_DURR_el$imps,
                     IURR_la = imp_IURR_la$imps,
                     IURR_el = imp_IURR_el$imps,
                     bridge  = imp_bridge$imps,
                     blasso  = imp_blasso$imps,
                     MI_PCA  = imp_PCA$mids,
                     MI_CART = imp_CART$imps,
                     MI_RF   = imp_RANF$imps,
                     MI_OP   = imp_MICE_OP$mids) # I need this for convergence
  
  ## Analyse --------------------------------------------------------------- ##
  # For each imp method, analyse all datasets based on model defined in init.R
  
  ## MLE estaimte of mean, variance, covariance
  
  # Multiple imputed dataset
  
  sem_fits <- lapply(list(DURR_la = imp_DURR_la$dats,
                          DURR_el = imp_DURR_el$dats,
                          IURR_la = imp_IURR_la$dats,
                          IURR_el = imp_IURR_el$dats,
                          bridge  = imp_bridge$dats,
                          blasso  = imp_blasso$dats,
                          MI_PCA  = imp_PCA$dats,
                          MI_CART = imp_CART$dats,
                          MI_RF   = imp_RANF$dats,
                          MI_OP   = imp_MICE_OP$dats), 
                     fit_sat_model)
  
  # Single dataset
  sem_sndt <- lapply(list(missFor = imp_missFor$dats,          
                          GS      = Xy,                        
                          CC      = Xy_mis[rowSums(!O) == 0, ]), 
                     sem, model = parms$lav_model, likelihood = "wishart")
  
  sem_sndt_output <- capture.output(sem_sndt)
  
  ## LM model
  
  # Multiple datasets
  
  lm_fits <- lapply(list(DURR_la = imp_DURR_la$dats,
                         DURR_el = imp_DURR_el$dats,
                         IURR    = imp_IURR_la$dats,
                         IURR_el = imp_IURR_el$dats,
                         bridge  = imp_bridge$dats,
                         blasso  = imp_blasso$dats,
                         MI_PCA  = imp_PCA$dats,
                         MI_CART = imp_CART$dats,
                         MI_RF   = imp_RANF$dats,
                         MI_OP   = imp_MICE_OP$dats), 
                    fit_lm_models, vrbs = parms$lm_model)
  
  # Single dataset
  
  lm_sndt <- fit_lm_models(list(missFor = imp_missFor$dats, 
                                GS      = Xy, 
                                CC      = Xy_mis[rowSums(!O) == 0, ]),
                           vrbs = parms$lm_model)
  
  ## Pooling --------------------------------------------------------------- ##
  # For each imp method, pool estimates across the m datasets
  
  # MLE mean, var, cov
  
  sem_pool_EST <- sapply(sem_fits[lapply(sem_fits, length) != 0], 
                         sem_pool_EST_f)
  sem_pool_CI  <- sapply(sem_fits[lapply(sem_fits, length) != 0], 
                         sem_pool_CI_f)
  
  # append single imputations, and GS and CC results
  sem_pool_EST <- cbind(sem_pool_EST,
                        sem_EST(sem_sndt))
  
  sem_pool_CI <- cbind(sem_pool_CI,
                       sem_CI(sem_sndt))
  
  # LM models
  lm_pool_est <- sapply(lm_fits[lapply(lm_fits, length) != 0], lm_pool_EST_f)
  lm_pool_CI <- sapply(lm_fits[lapply(lm_fits, length) != 0], lm_pool_CI_f)
  
  # append single imputations, and GS and CC results
  lm_pool_EST <- cbind(lm_pool_est,
                       lm_EST(lm_sndt))
  
  lm_pool_CI <- cbind(lm_pool_CI,
                      lm_CI(lm_sndt))
  
  ## Times ----------------------------------------------------------------- ##
  # aggregate imputation times
  
  imp_time <- list(DURR_la = imp_DURR_la$time,
                   DURR_el = imp_DURR_el$time,
                   IURR    = imp_IURR_la$time,
                   IURR_el = imp_IURR_el$time,
                   bridge  = imp_bridge$time,
                   blasso  = imp_blasso$time,
                   MI_PCA  = imp_PCA$time,
                   MI_CART = imp_CART$time,
                   MI_RF   = imp_RANF$time,
                   MI_OP   = imp_MICE_OP$time)
  imp_time <- do.call(cbind, imp_time)[1,]
  
  ## Store output ---------------------------------------------------------- ##
  output <- list(cond         = cond,
                 dat_full     = Xy,
                 dat_miss     = Xy_mis,
                 sem_EST      = sem_pool_EST,
                 sem_CI       = sem_pool_CI,
                 lm_EST       = lm_pool_EST,
                 lm_CI        = lm_pool_CI,
                 miss_descrps = miss_descrps,
                 run_time_min = imp_time,
                 imp_values   = imp_values)
  return(output)
}
