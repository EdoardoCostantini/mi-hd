### Title:    Convergence Checks: Initialization Script for Exp 2
### Project:  Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-08-11
### Notes:    The goal is to find out how many iterations should be used
###           in the full study for each method.

# Store these objects
parms$store <- c(cond     = TRUE,
                 dat_full = FALSE,
                 dat_miss = FALSE,
                 sem_EST  = FALSE,
                 sem_CI   = FALSE,
                 CFA_EST  = FALSE,
                 CFA_CI   = FALSE,
                 semS_EST = FALSE,
                 semS_CI  = FALSE,
                 fmi      = TRUE,
                 miss_des = FALSE,
                 time     = TRUE,
                 time_prep = TRUE,
                 imps     = TRUE)

# Decide which methods you car about for the convergence check
parms$meth_sel <- list(DURR_all = FALSE,   # version w/o SI
                       DURR_si  = FALSE,  # version w/o SI
                       IURR_all = TRUE,   # version w/o SI
                       IURR_si  = TRUE,  # version w/ SI
                       blasso   = FALSE,
                       bridge   = FALSE,
                       MI_PCA   = FALSE,
                       MI_CART  = FALSE,
                       MI_RF    = FALSE,
                       MI_OP    = FALSE,
                       missFor  = TRUE,
                       mean     = TRUE,
                       CC       = TRUE,
                       GS       = TRUE)

parms$methods <- names(parms$meth_sel)[which(parms$meth_sel==TRUE)]


# Always keep missForest to a minimal computation time
parms$missFor_maxiter <- 2 # 20
parms$missFor_ntree   <- 10 # 100

# Decide which paramters suit best the convergence check
parms$dt_rep     <- 5 # 5
parms$chains     <- 4 # 4
parms$iters      <- 250 # 250 if you do not converge before this, the method is not feasable 
parms$burnin_imp <- 0 # no need, I want to see the first part as well
parms$ndt        <- parms$chains # same chians as dataset in mice like imp
parms$thin       <- (parms$iters - parms$burnin_imp)/parms$ndt
parms$pos_dt     <- (parms$burnin_imp+1):parms$iters # candidate datasets (after convergence)
parms$keep_dt    <- parms$pos_dt[seq(1, length(parms$pos_dt), parms$thin)] # keep 1 dataset every thin

# For blasso
parms$chains_bl     <- parms$chains
parms$iters_bl      <- 250 # 250
parms$burnin_imp_bl <- 0 # 0
parms$thin_bl       <- (parms$iters_bl - parms$burnin_imp_bl)/parms$ndt
parms$pos_dt_bl     <- (parms$burnin_imp_bl+1):parms$iters_bl # candidate datasets
parms$keep_dt_bl    <- parms$pos_dt_bl[seq(1, 
                                           length(parms$pos_dt_bl), 
                                           parms$thin_bl)]

# For mice-like algorithms
parms$mice_iters <- parms$iters # for this they are the same.
parms$mice_ndt   <- parms$ndt # 10 # number of imputed datasets to pool esitmaes from (10)

# Change name of outputs

# Output and Progres report related
parms$outDir <- "../output/"
parms$start_time <- format(Sys.time(), "%Y%m%d_%H%M")

parms$report_file_name <- paste0("exp",
                                 parms$exp, "_",
                                 "conv_",
                                 parms$start_time, 
                                 ".txt")
parms$results_file_name <- paste0("exp",
                                  parms$exp, "_",
                                  "conv_",
                                  parms$start_time,
                                  ".rds")
parms$description <- c("For the most challanging condition, data generation is repeated 
                        a few times, imputation models are used. Obtained imputations can 
                        then be used to check convergence")

# Fix condition to desired one for convergence check
conds <- conds[nrow(conds), ]
