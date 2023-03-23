# Title:    Imputing High Dimensional Data
# Author:   Edoardo Costantini
# Created:  2020-05-19
# Modified: 2023-03-22

# Functions ---------------------------------------------------------------

## Run one replication of the simulation:
doRep <- function(rp, conds, parms, debug = FALSE, verbose = TRUE) {
  ## For internals
  # rp = 1
  # debug = FALSE
  # verbose = FALSE
  # i = 1
  
  ## For 1 repetition, performs simulation across all conditions
  ## Setup the PRNG: (for setting seed with guaranteed independece)
  # form the rlecuyer package
  # rp is the replication index, in the mcapply function it the argument
  # that takes the elements of the vector specified for X
  # This creates names for the streams so that you can recover them
  
  ## START capture output
  if(verbose == TRUE){
    sink(paste0(parms$outDir, "console/", 
                parms$report_file_name, "_dataRep", rp, ".txt"))
  }
  
  ## Set seed
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
    names(rp_out) <- paste0("cond_", 
                            sapply(1:nrow(conds), 
                                   function(x) paste0(conds[x,], 
                                                      collapse = "_")
                                   )
                            )
    
# Experiment 1 ------------------------------------------------------------

    if(parms$exp == 1){
      for(i in 1 : nrow(conds)) {
        
        print(paste0("Exp: ", parms$exp, 
                     " > Rep: ", rp, 
                     " > Cond: ", i,
                     " - Time: ",
                     Sys.time()))
        
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
    }

# Experiment 2 ------------------------------------------------------------
    
    if(parms$exp == 2){
      for(i in 1 : nrow(conds)) {
        
        print(paste0("Exp: ", parms$exp, 
                     " > Rep: ", rp, 
                     " > Cond: ", i,
                     " - Time: ",
                     Sys.time()))
        
        if(debug == TRUE){
          rp_out[[i]] <- capture.output(tryCatch({
            withErrorTracing({runCell_lv(cond = conds[i, ],
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
              runCell_lv(cond = conds[i, ],
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
    }
    
# Experiment 4 ------------------------------------------------------------
    
    if(parms$exp == 4){
      for(i in 1 : nrow(conds)) {
        
        print(paste0("--- Exp: ", parms$exp, 
                     " > Rep: ", rp, 
                     " > Cond: ", i,
                     " - Time: ",
                     Sys.time(),
                     " ---"))
        
        if(debug == TRUE){
          rp_out[[i]] <- capture.output(tryCatch({
            withErrorTracing({runCell_evs(cond = conds[i, ],
                                          parms = parms,
                                          rep_status = rp,
                                          data_source = readRDS("../data/exp4_EVS2017_full.rds")$full)})
          }, error = function(e){
            e <<- e
            cat("ERROR: ", e$message, "\nin ")
            print(e$call)
          }))
        } else {
          rp_out[[i]] <- tryCatch(
            {
              # Try running simulation for condition i, repetition rp
              runCell_evs(cond = conds[i, ],
                          parms = parms,
                          rep_status = rp,
                          data_source = readRDS("../data/exp4_EVS2017_full.rds")$full)
            },
            error = function(report) {
              err <- paste0("Original Error: ", report)
              print(err)
              return(err)
            }
          )
        }
      }
    }

  ## END capture output
  if(verbose == TRUE){
    sink()
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
  # source("./exp1_init.R")
  # rep_status <- 1
  # cond <- conds[4, ]
  
  ## Data ------------------------------------------------------------------ ##
  # According to experiment set up, gen 1 fully-obs data dataset and
  # impose missing values;
  
  Xy <- simData_exp1(cond, parms)
  Xy_mis <- imposeMiss(Xy, parms, cond)
  Xy_mis <- cbind(Xy_mis[, parms$z_m_id], 
                  Xy_mis[, -which(colnames(Xy_mis) %in% parms$z_m_id)])
  
  O <- !is.na(Xy_mis) # matrix index of observed values
  miss_descrps <- colMeans(!O[, parms$z_m_id]) 

  ## Imputation ------------------------------------------------------------ ##
  # Impute m times the data w/ missing values w/ different methods

  # Impute according to DURR method
  
    # Lasso
  imp_DURR_la <- impute_DURR(Z = Xy_mis,
                             O = as.data.frame(O),
                             reg_type = "lasso",
                             cond = cond,
                             perform = parms$meth_sel$DURR_la,
                             parms = parms)

  # Impute according to IURR method
    # Lasso
  imp_IURR_la <- impute_IURR(Z = Xy_mis,
                             O = as.data.frame(O),
                             reg_type = "lasso",
                             cond = cond,
                             perform = parms$meth_sel$IURR_la,
                             parms = parms)

  # Impute according to Hans Blasso method
  imp_blasso <- impute_BLAS_hans(Z = Xy_mis,
                                 parms = parms,
                                 perform = parms$meth_sel$blasso)

  # Impute according to van Buuren Ridge
  imp_bridge <- impute_BRIDGE(Z = Xy_mis, 
                              O = as.data.frame(O),
                              ridge_p = cond$ridge,
                              parms = parms,
                              perform = parms$meth_sel$bridge)

  # Impute according to Howard Et Al 2015 PCA appraoch
  imp_PCA <- impute_PCA(Z = Xy_mis, O = O, cond = cond, parms = parms)

  # MICE-CART traditional
  imp_CART <- impute_CART(Z = Xy_mis,
                          O = as.data.frame(O),
                          perform = parms$meth_sel$MI_CART,
                          parms = parms)
  
  # MICE-RF
  imp_RANF <- impute_RANF(Z = Xy_mis,
                          O = as.data.frame(O),
                          perform = parms$meth_sel$MI_RF,
                          parms = parms)

  # MICE w/ step-wise forward regression
  imp_stepFor <- impute_IVEware(
    Z = Xy_mis,
    minR2 = cond$minR2,
    rep_status = rep_status,
    parms = parms
  )

  # MICE w/ quickpred
  imp_qp <- impute_MICE_qp(Z = Xy_mis,
                           perform = parms$meth_sel$MI_qp,
                           ridge = 1e-5, eps = 1e-04, threshold = .999,
                           parms = parms)

  # MICE w/ analysis model only
  imp_MICE_am <- impute_MICE_cp(Z = Xy_mis,
                                preds = parms$z_m_id,
                                perform = parms$meth_sel$MI_am,
                                ridge = 1e-5, eps = 1e-04, threshold = .999,
                                parms = parms)

  # MICE w/ true model
  imp_MICE_OP <- impute_MICE_OP(Z = Xy_mis,
                                O = O,
                                cond = cond,
                                perform = parms$meth_sel$MI_OP,
                                parms = parms)

  # missForest
  imp_missFor <- impute_missFor(Z = Xy_mis, parms = parms)
  
  # mean imputation
  Xy_mean <- Xy_mis
  for (j in 1:parms$zm_n) {
    ry <- !is.na(Xy_mean[, parms$z_m_id[j]])
    mu_zj <- mean(Xy_mean[, parms$z_m_id[j]], na.rm = TRUE)
    Xy_mean[!ry, parms$z_m_id[j]] <- mu_zj
  }
  
  ## Convergence ----------------------------------------------------------- ##
  
  imp_values <- list(DURR_la = imp_DURR_la$imps,
                     IURR_la = imp_IURR_la$imps,
                     bridge  = imp_bridge$imps,
                     blasso  = imp_blasso$imps,
                     MI_PCA  = imp_PCA$mids,
                     MI_CART = imp_CART$imps,
                     MI_RF   = imp_RANF$imps,
                     stepFor = imp_stepFor$imps,
                     MI_qp   = imp_qp$imps,
                     MI_am   = imp_MICE_am$imps,
                     MI_OP   = imp_MICE_OP$mids) # I need this for convergence
  
  ## Analyse --------------------------------------------------------------- ##
  # For each imp method, analyse all datasets based on model defined in init.R
  
  ## MLE estimate of mean, variance, covariance
  
  # Multiple imputed dataset
  
  sem_fits <- lapply(list(DURR_la = imp_DURR_la$dats,
                          IURR_la = imp_IURR_la$dats,
                          bridge  = imp_bridge$dats,
                          blasso  = imp_blasso$dats,
                          MI_PCA  = imp_PCA$dats,
                          MI_CART = imp_CART$dats,
                          MI_RF   = imp_RANF$dats,
                          stepFor = imp_stepFor$dats,
                          MI_qp   = imp_qp$dats,
                          MI_am   = imp_MICE_am$dats,
                          MI_OP   = imp_MICE_OP$dats), 
                     fit_sat_model)
  
  # Single dataset
  sem_sndt <- lapply(list(missFor = imp_missFor$dats,          
                          mean    = Xy_mean,
                          CC      = Xy_mis[rowSums(!O) == 0, ],
                          GS      = Xy), 
                     sem, model = parms$lav_model, likelihood = "wishart")
  
  ## LM model
  
  # Multiple datasets
  
  lm_fits <- lapply(list(DURR_la = imp_DURR_la$dats,
                         IURR_la = imp_IURR_la$dats,
                         bridge  = imp_bridge$dats,
                         blasso  = imp_blasso$dats,
                         MI_PCA  = imp_PCA$dats,
                         MI_CART = imp_CART$dats,
                         MI_RF   = imp_RANF$dats,
                         stepFor = imp_stepFor$dats,
                         MI_qp   = imp_qp$dats,
                         MI_am   = imp_MICE_am$dats,
                         MI_OP   = imp_MICE_OP$dats), 
                    fit_lm_models, mod = parms$lm_model)
  
  # Single dataset
  
  lm_sndt <- fit_lm_models(list(missFor = imp_missFor$dats,          
                                mean    = Xy_mean,
                                CC      = Xy_mis[rowSums(!O) == 0, ],
                                GS      = Xy),
                           mod = parms$lm_model)
  
  ## Pooling --------------------------------------------------------------- ##
  # For each imp method, pool estimates across the m datasets
  
  # MLE mean, var, cov
  
  sem_pool_MI_EST <- sapply(sem_fits[lapply(sem_fits, length) != 0], 
                            sem_pool_EST_f)
  sem_pool_MI_CI  <- sapply(sem_fits[lapply(sem_fits, length) != 0], 
                            sem_pool_CI_f)
  sem_fmi <- sapply(sem_fits[lapply(sem_fits, length) != 0], 
                    .fmi_compute, form = "finite")
  
  # append single imputations, and GS and CC results
  sem_gather_EST <- cbind(sem_pool_MI_EST,
                          sem_EST(sem_sndt))
  
  sem_gather_CI <- cbind(sem_pool_MI_CI,
                         sem_CI(sem_sndt))
  
  # LM models
  lm_pool_est <- sapply(lm_fits[lapply(lm_fits, length) != 0], lm_pool_EST_f)
  lm_pool_CI <- sapply(lm_fits[lapply(lm_fits, length) != 0], lm_pool_CI_f)
  lm_fmi <- sapply(lm_fits[lapply(lm_fits, length) != 0], 
                    .fmi_compute, form = "finite")
  
  # append single imputations, and GS and CC results
  lm_pool_EST <- cbind(lm_pool_est,
                       lm_EST(lm_sndt))
  
  lm_pool_CI <- cbind(lm_pool_CI,
                      lm_CI(lm_sndt))
  
  
  ## Times ----------------------------------------------------------------- ##
  # aggregate imputation times
  
  imp_time <- list(DURR_la = imp_DURR_la$time,
                   IURR_la = imp_IURR_la$time,
                   bridge  = imp_bridge$time,
                   blasso  = imp_blasso$time,
                   MI_PCA  = imp_PCA$time,
                   MI_CART = imp_CART$time,
                   MI_RF   = imp_RANF$time,
                   stepFor = imp_stepFor$time,
                   MI_qp   = imp_qp$time,
                   MI_am   = imp_MICE_am$time,
                   MI_OP   = imp_MICE_OP$time)
  imp_time <- do.call(cbind, imp_time)[1,]
  
  ## Store output ---------------------------------------------------------- ##
  output <- list(cond         = cond,
                 dat_full     = Xy,
                 dat_miss     = Xy_mis,
                 sem_EST      = sem_gather_EST,
                 sem_CI       = sem_gather_CI,
                 lm_EST       = lm_pool_EST,
                 lm_CI        = lm_pool_CI,
                 fmi          = list(sems = sem_fmi,
                                     lm   = lm_fmi),
                 miss_descrps = miss_descrps,
                 run_time_min = imp_time,
                 imp_values   = imp_values)[parms$store]
  
  ## Return Output of function
  return(output)
}

# Experiment 2 ------------------------------------------------------------

runCell_lv <- function(cond, parms, rep_status) {
  ## Description
  # Given 1 condition, Generates 1 dataset and performs imutations according to
  # selected methods
  ## For internals
  # set.seed(1234)
  # cond <- conds[2, ]
  
  ## ----------------------------------------------------------------------- ##
  
  # ------------------- #
    ## Data ####
  # ------------------- #
  
  simData_list <- simData_lv(parms, cond)
    Xy <- simData_list$dat
  Xy_mis <- imposeMiss_lv(simData_list, parms, cond)
  
  O <- !is.na(Xy_mis) # matrix index of observed values
  miss_descrps <- colMeans(!O[, 1:parms$zm_n]) 
  
  ## ----------------------------------------------------------------------- ##
  
  # ------------------- #
    ## Imputation ####
  # ------------------- #
  
  # Impute m times the data w/ missing values w/ different methods
  
  ## ----------------------------------------------------------------------- ##
  # Impute according to DURR method
  # Lasso
  imp_DURR_la <- impute_DURR(Z = Xy_mis,
                             O = as.data.frame(O),
                             reg_type = "lasso",
                             cond = cond,
                             perform = parms$meth_sel$DURR_la,
                             parms = parms)

  ## ----------------------------------------------------------------------- ##
  # Impute according to IURR method
  # Lasso
  imp_IURR_la <- impute_IURR(Z = Xy_mis,
                             O = as.data.frame(O),
                             reg_type = "lasso",
                             cond = cond,
                             perform = parms$meth_sel$IURR_la,
                             parms = parms)
  
  ## ----------------------------------------------------------------------- ##
  # Impute according to Hans Blasso method
  imp_blasso <- impute_BLAS_hans(Z = Xy_mis, 
                                 # O = as.data.frame(O),
                                 parms = parms,
                                 perform = parms$meth_sel$blasso)

  ## ----------------------------------------------------------------------- ##
  # Impute according to van Buuren Ridge
  imp_bridge <- impute_BRIDGE(Z = Xy_mis, 
                              O = as.data.frame(O),
                              ridge_p = cond$ridge,
                              parms = parms,
                              perform = parms$meth_sel$bridge)

  ## ----------------------------------------------------------------------- ##
  # Impute according to Howard Et Al 2015 PCA appraoch
  imp_PCA <- impute_PCA(Z = Xy_mis, O = O, cond = cond, parms = parms)

  ## ----------------------------------------------------------------------- ##
  # MICE-CART traditional
  imp_CART <- impute_CART(Z = Xy_mis,
                          O = as.data.frame(O),
                          cond = cond,
                          perform = parms$meth_sel$MI_CART,
                          parms = parms)

  ## ----------------------------------------------------------------------- ##
  # MICE-RF
  imp_RANF <- impute_RANF(Z = Xy_mis,
                          O = as.data.frame(O),
                          cond = cond,
                          perform = parms$meth_sel$MI_RF,
                          parms = parms)

  ## ----------------------------------------------------------------------- ##
  # MICE w/ true model
  imp_MICE_OP <- impute_MICE_OP(Z = Xy_mis,
                                O = O,
                                cond = cond,
                                perform = parms$meth_sel$MI_OP,
                                parms = parms)
  
  ## ----------------------------------------------------------------------- ##
  # missForest
  imp_missFor <- impute_missFor(Z = Xy_mis, parms = parms)
  
  ## ----------------------------------------------------------------------- ##
  # mean imputation
  Xy_mean <- Xy_mis
  for (j in 1:parms$zm_n) {
    ry <- !is.na(Xy_mean[, parms$z_m_id[j]])
    mu_zj <- mean(Xy_mean[, parms$z_m_id[j]], na.rm = TRUE)
    Xy_mean[!ry, parms$z_m_id[j]] <- mu_zj
  }
  
  ## ----------------------------------------------------------------------- ##
  
  # ------------------- #
    ## Convergence ####
  # ------------------- #
  
  imp_values <- list(DURR_la = imp_DURR_la$imps,
                     IURR_la = imp_IURR_la$imps,
                     bridge  = imp_bridge$imps,
                     blasso  = imp_blasso$imps,
                     MI_PCA  = imp_PCA$mids,
                     MI_CART = imp_CART$imps,
                     MI_RF   = imp_RANF$imps,
                     MI_OP   = imp_MICE_OP$mids) # I need this for convergence
  
  ## ----------------------------------------------------------------------- ##
  
  # ------------------- #
  ## Analyse and pool ####
  # ------------------- #
  
  # For each imp method, analyse all datasets based on model defined in init.R
  
  ## Create Scored data (needed to define analysis model)
  SC_dt_sn <- lapply(list(missFor = imp_missFor$dats,          
                          mean    = Xy_mean,
                          CC      = Xy_mis[rowSums(!O) == 0, ],
                          GS      = Xy),
                     scorify, cond = cond, parms = parms)
  
  SC_dt_mi <-  lapply(list(DURR_la = imp_DURR_la$dats,
                           IURR_la = imp_IURR_la$dats,
                           bridge  = imp_bridge$dats,
                           blasso  = imp_blasso$dats,
                           MI_PCA  = imp_PCA$dats,
                           MI_CART = imp_CART$dats,
                           MI_RF   = imp_RANF$dats,
                           MI_OP   = imp_MICE_OP$dats), 
                      function(x){
                        lapply(x, scorify, cond = cond, parms = parms)
                      }
  )
  
  ## Define Analysis models
  SAT_mod_raw <- SAT_mod_write(parms$z_m_id) # raw data
  CFA_mod_raw <- CFA_mod_wirte(Xy, 2, parms) # raw data
  SAT_mod_sco <- SAT_mod_write(colnames(SC_dt_sn$GS)[1:parms$sc_n]) # score data
  lm_formula <- "sc1 ~ -1 + sc2 + sc3"
  
  ## ----------------------------------------------------------------------- ##
  
  # ------------------- #
    ## (1) SEM raw ####
  # ------------------- #
    
  # Fit models
  semR_fit_mi <- lapply(list(DURR_la = imp_DURR_la$dats,
                             IURR_la = imp_IURR_la$dats,
                             bridge  = imp_bridge$dats,
                             blasso  = imp_blasso$dats,
                             MI_PCA  = imp_PCA$dats,
                             MI_CART = imp_CART$dats,
                             MI_RF   = imp_RANF$dats,
                             MI_OP   = imp_MICE_OP$dats), 
                        fit_sem, model = SAT_mod_raw)
  
  semR_fit_sn <- fit_sem(list(missFor = imp_missFor$dats,          
                              mean    = Xy_mean,
                              CC      = Xy_mis[rowSums(!O) == 0, ],
                              GS      = Xy),
                         model = SAT_mod_raw)
  
  # Pool paramters
  semR_est <- sapply(semR_fit_mi[lapply(semR_fit_mi, length) != 0], 
                     sem_pool_EST_f)
  semR_ci  <- sapply(semR_fit_mi[lapply(semR_fit_mi, length) != 0], 
                     sem_pool_CI_f)
  semR_fmi <- sapply(semR_fit_mi[lapply(semR_fit_mi, length) != 0], 
                    .fmi_compute)
  
  # Prep final (append single imputations, GS and CC results)
  semR_est_all <- cbind(semR_est,
                        sem_EST(semR_fit_sn))
  
  semR_ci_all <- cbind(semR_ci,
                       sem_CI(semR_fit_sn))
  
  ## ----------------------------------------------------------------------- ##
  
  # ------------------- #
  ## (2) CFA ####
  # ------------------- #
    
  # Fit models
  CFA_fit_mi <- lapply(list(DURR_la = imp_DURR_la$dats,
                            IURR_la = imp_IURR_la$dats,
                            bridge  = imp_bridge$dats,
                            blasso  = imp_blasso$dats,
                            MI_PCA  = imp_PCA$dats,
                            MI_CART = imp_CART$dats,
                            MI_RF   = imp_RANF$dats,
                            MI_OP   = imp_MICE_OP$dats), 
                       fit_sem, model = CFA_mod_raw, std.lv = TRUE)
  CFA_fit_sn <- fit_sem(list(missFor = imp_missFor$dats,          
                             mean    = Xy_mean,
                             CC      = Xy_mis[rowSums(!O) == 0, ],
                             GS      = Xy),
                        model = CFA_mod_raw,
                        std.lv = TRUE)

  # Pool parameters
  CFA_est_mi <- sapply(CFA_fit_mi[lapply(CFA_fit_mi, length) != 0], 
                            sem_pool_EST_f)
  CFA_ci_mi  <- sapply(CFA_fit_mi[lapply(CFA_fit_mi, length) != 0], 
                            sem_pool_CI_f)
  CFA_fmi <- sapply(CFA_fit_mi[lapply(CFA_fit_mi, length) != 0], 
                    .fmi_compute)
  
  # Prep final (append single imputations, GS and CC results)
  CFA_est_all <- cbind(CFA_est_mi,
                       sem_EST(CFA_fit_sn))
  CFA_ci_all <- cbind(CFA_ci_mi,
                      sem_CI(CFA_fit_sn))
  
  ## ----------------------------------------------------------------------- ##
  
  # ------------------- #
  ## (3) SEM on scored data ####
  # ------------------- #
    
  # Fit models
  semS_fit_mi <- lapply(SC_dt_mi,
                          fit_sem, model = SAT_mod_sco)
  semS_fit_sn <- fit_sem(SC_dt_sn,
                           model = SAT_mod_sco)
  
  # Pool parameters (and print FMI)
  semS_est <- sapply(semS_fit_mi[lapply(semS_fit_mi, length) != 0], 
                       sem_pool_EST_f)[-c(3, 6),]
  semS_ci  <- sapply(semS_fit_mi[lapply(semS_fit_mi, length) != 0], 
                       sem_pool_CI_f)[-c(3, 6, 12, 15),]
  semS_fmi <- sapply(semS_fit_mi[lapply(semS_fit_mi, length) != 0], 
                     .fmi_compute)
  
  # Prep final (append single imputations, GS and CC results)
  semS_est_all <- cbind(semS_est,
                        sem_EST(semS_fit_sn)[-c(3, 6),])
  
  semS_ci_all <- cbind(semS_ci,
                       sem_CI(semS_fit_sn)[-c(3, 6, 12, 15),])
  
  
  ## ----------------------------------------------------------------------- ##
  
  # ------------------- #
  ## (4) LM model ####
  # ------------------- #
  
  # Fit models
  lm_mi <- lapply(SC_dt_mi, 
                         fit_lm, model = lm_formula)
  lm_sn <- lapply(SC_dt_sn, lm, formula = lm_formula)
  
  # Pool paramters
  lm_est <- sapply(lm_mi[lapply(lm_mi, length) != 0], 
                            lm_pool_EST_f)
  lm_ci  <- sapply(lm_mi[lapply(lm_mi, length) != 0], 
                            lm_pool_CI_f)
  lm_fmi <- sapply(lm_mi[lapply(lm_mi, length) != 0], 
                     .fmi_compute)
  
  # Prep final (append single imputations, GS and CC results)
  lm_est_all <- cbind(lm_est,
                      lm_EST(lm_sn))
  
  lm_ci_all <- cbind(lm_ci,
                      lm_CI(lm_sn))
  
  ## ----------------------------------------------------------------------- ##
  
  # ------------------- #
    ## Times ####
  # ------------------- #
    
  # aggregate imputation times
  
  imp_time <- list(DURR_la = imp_DURR_la$time,
                   IURR_la = imp_IURR_la$time,
                   bridge  = imp_bridge$time,
                   blasso  = imp_blasso$time,
                   MI_PCA  = imp_PCA$time,
                   MI_CART = imp_CART$time,
                   MI_RF   = imp_RANF$time,
                   MI_OP   = imp_MICE_OP$time)
  imp_time <- do.call(cbind, imp_time)[1,]
  
  ## ----------------------------------------------------------------------- ##
  
  # ------------------- #
    ## Store Output ####
  # ------------------- #
  
  output <- list(cond         = cond,
                 dat_full     = Xy,
                 dat_miss     = Xy_mis,
                 # SEM raw
                 semR_EST     = semR_est_all,
                 semR_CI      = semR_ci_all,
                 # CFA raw
                 CFA_EST      = CFA_est_all,
                 CFA_CI       = CFA_ci_all,
                 # SEM Scored
                 semS_EST     = semS_est_all,
                 semS_CI      = semS_ci_all,
                 # LM Scored
                 lm_EST       = lm_est_all,
                 lm_CI        = lm_ci_all,
                 # Other
                 fmi          = list(semR = semR_fmi,
                                     CFA  = CFA_fmi,
                                     semS = semS_fmi,
                                     lm   = lm_fmi),
                 miss_descrps = miss_descrps,
                 run_time_min = imp_time,
                 imp_values   = imp_values)[parms$store]
  
  return(output)
  
}

# Experiment 4 ------------------------------------------------------------

runCell_evs <- function(cond, parms, rep_status, data_source) {
  ## Description
  # Given 1 condition, Generates 1 dataset and performs imputations 
  # according to selected methods
  ## For internals
  # set.seed(1234)
  # cond <- conds[1, ]
  # data_source <- readRDS("../data/exp4_EVS2017_full.rds")$full
  
  ## Data ------------------------------------------------------------------ ##
  # Gen one fully-obs data
  Xy <- data_source[sample(1:nrow(data_source),
                           cond$n,
                           replace = TRUE), ]
  
  # Impose missing values
  Xy_mis <- imposeMiss_evs(Xy, parms, cond)
  
  # Missing data
  miss_descrps <- colMeans(is.na(Xy_mis)[, parms$z_m_id])
  
  ## Imputation ------------------------------------------------------------ ##
  # Impute m times the data w/ missing values w/ different methods
  
  # Impute according to DURR method
  # Lasso
  imp_DURR_la <- impute_DURR(Z = Xy_mis,
                             O = data.frame(!is.na(Xy_mis)),
                             reg_type = "lasso",
                             cond = cond,
                             perform = parms$meth_sel$DURR_la,
                             parms = parms)
  
  # Impute according to IURR method
  # Lasso
  imp_IURR_la <- impute_IURR(Z = Xy_mis,
                             O = data.frame(!is.na(Xy_mis)),
                             reg_type = "lasso",
                             cond = cond,
                             perform = parms$meth_sel$IURR_la,
                             parms = parms)
  
  # Impute according to Hans Blasso method
  imp_blasso <- impute_BLAS_hans(Z = Xy_mis,
                                 parms = parms,
                                 perform = parms$meth_sel$blasso)
  
  # Impute according to van Buuren Ridge
  imp_bridge <- impute_BRIDGE(Z = Xy_mis,
                              O = data.frame(!is.na(Xy_mis)),
                              ridge_p = cond$ridge,
                              parms = parms,
                              perform = parms$meth_sel$bridge)
  
  # Impute according to Howard Et Al 2015 PCA appraoch
  imp_PCA <- impute_PCA(Z     = Xy_mis,
                        O     = data.frame(!is.na(Xy_mis)),
                        cond  = cond,
                        parms = parms)
  
  # MICE-CART
  imp_CART <- impute_CART(Z = Xy_mis,
                          O = data.frame(!is.na(Xy_mis)),
                          cond = cond,
                          perform = parms$meth_sel$MI_CART,
                          parms = parms)
  
  # MICE-RF
  imp_RANF <- impute_RANF(Z = Xy_mis,
                          O = data.frame(!is.na(Xy_mis)),
                          cond = cond,
                          perform = parms$meth_sel$MI_RF,
                          parms = parms)


  # MICE w/ step-wise forward regression
  imp_stepFor <- impute_IVEware(
    Z = Xy_mis,
    minR2 = cond$minR2,
    rep_status = rep_status,
    parms = parms
  )


  # MICE w/ quickpred
  imp_qp <- impute_MICE_qp(Z = Xy_mis,
                           perform = parms$meth_sel$MI_qp,
                           ridge = 1e-5, eps = 1e-04, threshold = .999,
                           # mice defaults!
                           parms = parms)

  # MICE w/ analysis model only
  imp_MICE_am <- impute_MICE_cp(Z = Xy_mis,
                                preds = parms$am_vars,
                                ridge = 1e-5, eps = 1e-04, threshold = .999,
                                # mice defaults!
                                perform = parms$meth_sel$MI_am,
                                parms = parms)

  # MICE w/ true model
  imp_MICE_OP <- impute_MICE_OP(Z = Xy_mis,
                                O = data.frame(!is.na(Xy_mis)),
                                cond = cond,
                                perform = parms$meth_sel$MI_OP,
                                parms = parms)
  
  # missForest
  imp_missFor <- impute_missFor(Z = Xy_mis, parms = parms)
  
  # mean imputation
  Xy_mean <- Xy_mis
  for (j in 1:parms$zm_n) {
    ry <- !is.na(Xy_mean[, parms$z_m_id[j]])
    mu_zj <- mean(Xy_mean[, parms$z_m_id[j]], na.rm = TRUE)
    Xy_mean[!ry, parms$z_m_id[j]] <- mu_zj
  }
  
## Convergence ------------------------------------------------------------- ##

  imp_values <- list(DURR_la    = imp_DURR_la$imps,
                     IURR_la    = imp_IURR_la$imps,
                     blasso     = imp_blasso$imps,
                     bridge     = imp_bridge$imps,
                     MI_PCA     = imp_PCA$mids,
                     MI_CART    = imp_CART$imps,
                     MI_RF      = imp_RANF$imps,
                     stepFor    = imp_stepFor$imps,
                     MI_qp      = imp_qp$imps,
                     MI_am      = imp_MICE_am$imps,
                     MI_OP      = imp_MICE_OP$mids)

## Analyse ----------------------------------------------------------------- ##
  # For each imp method, analyse all datasets based on model defined in init.R
  
  # Linear Models
  # Data
  mi_data <- list(DURR_la    = imp_DURR_la$dats,
                  IURR_la    = imp_IURR_la$dats,
                  blasso     = imp_blasso$dats,
                  bridge     = imp_bridge$dats,
                  MI_PCA     = imp_PCA$dats,
                  MI_CART    = imp_CART$dats,
                  MI_RF      = imp_RANF$dats,
                  stepFor    = imp_stepFor$dats,
                  MI_qp      = imp_qp$dats,
                  MI_am      = imp_MICE_am$dats,
                  MI_OP      = imp_MICE_OP$dats)
  
  si_data <- list(missFor = imp_missFor$dats,
                  mean    = Xy_mean,
                  CC      = na.omit(Xy_mis),
                  GS      = Xy)                    
  
## Fit --------------------------------------------------------------------- ##
  
  # model 1
  m1_mi <- lapply(mi_data, exp4_fit_mod1)
  m1_sn <- exp4_fit_mod1(si_data)

  # Model 2
  m2_mi <- lapply(mi_data, exp4_fit_mod2)
  m2_sn <- exp4_fit_mod2(si_data)

## Pool MI paramters ------------------------------------------------------- ##
  
  # Model 1
  # Est
  m1_est.mi <- sapply(m1_mi[lapply(m1_mi, length) != 0], 
                   lm_pool_EST_f)
  m1_est.sn <- sapply(m1_sn, function(x) coef(x)[parms$m1_par])
  m1_est <- data.frame(cbind(m1_est.mi, m1_est.sn))
  
  # Confint
  m1_ci.mi  <- sapply(m1_mi[lapply(m1_mi, length) != 0], 
                   lm_pool_CI_f)
  m1_ci.sn  <- sapply(m1_sn, function(x){
    CI.temp <- confint(x)
    return(CI.temp[parms$m1_par, ])
  })
  m1_ci <- data.frame( cbind(m1_ci.mi, m1_ci.sn) )
  
  # FMI
  m1_fmi <- try(
    apply(sapply(m1_mi[lapply(m1_mi, length) != 0],
                 .fmi_compute), 2, mean, na.rm = TRUE),
    silent = TRUE)
  
  if(class(m1_fmi) == "try-error"){
    # If for some reason, a bridge value in the crossvalidation
    # preprocessing procedure does not do the trick, it is
    # automatically excluded by this if statement
    # I ideally this should be done in the result function not here!
    m1_fmi <- NA
  }
  
  # Model 2
  # Est
  m2_est.mi <- sapply(m2_mi[lapply(m2_mi, length) != 0], 
                      lm_pool_EST_f)
  m2_est.sn <- sapply(m2_sn, function(x) coef(x)[parms$m2_par])
  m2_est <- data.frame(cbind(m2_est.mi, m2_est.sn))
  
  # Confint
  m2_ci.mi  <- sapply(m2_mi[lapply(m2_mi, length) != 0], 
                      lm_pool_CI_f)
  m2_ci.sn  <- sapply(m2_sn, function(x){
    CI.temp <- confint(x)
    return(CI.temp[parms$m2_par, ])
  })
  m2_ci <- data.frame( cbind(m2_ci.mi, m2_ci.sn) )
  
  # FMI
  m2_fmi <- try(
    apply(sapply(m2_mi[lapply(m2_mi, length) != 0],
                 .fmi_compute), 2, mean, na.rm = TRUE),
    silent = TRUE)
  
  if(class(m2_fmi) == "try-error"){
    # If for some reason, a bridge value in the crossvalidation
    # preprocessing procedure does not do the trick, it is
    # automatically excluded by this if statement
    # I ideally this should be done in the result function not here!
    m2_fmi <- NA
  }

## Times ------------------------------------------------------------------- ##
  
  # aggregate imputation times
  imp_time <- list(DURR_la    = imp_DURR_la$time,
                   IURR_la    = imp_IURR_la$time,
                   blasso     = imp_blasso$time,
                   bridge     = imp_bridge$time,
                   MI_PCA     = imp_PCA$time,
                   MI_CART    = imp_CART$time,
                   MI_RF      = imp_RANF$time,
                   stepFor    = imp_stepFor$time,
                   MI_qp      = imp_qp$time,
                   MI_am      = imp_MICE_am$time,
                   MI_OP      = imp_MICE_OP$time)
  imp_time <- do.call(cbind, imp_time)[1,]
  
## Store output ------------------------------------------------------------ ##
  
  output <- list(cond         = cond,
                 # Datasets
                 dat_full     = Xy,
                 dat_miss     = Xy_mis,
                 # Linear model 1
                 m1_EST       = m1_est,
                 m1_CI        = m1_ci,
                 # Linear model 2
                 m2_EST       = m2_est,
                 m2_CI        = m2_ci,
                 # Other
                 fmi          = list(m1 = m1_fmi,
                                     m2 = m2_fmi),
                 miss_descrps = miss_descrps,
                 run_time_min = imp_time,
                 imp_values   = imp_values)[parms$store]
  return(output)
}