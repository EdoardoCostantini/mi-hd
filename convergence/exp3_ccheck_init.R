### Title:    Convergence Checks: Initialization Script for Exp 2
### Project:  Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-08-11
### Notes:    The goal is to find out how many iterations should be used
###           in the full study for each method.

# Get uninteresting parameters from the init_exp1.R script
source("./exp3_init.R")

# What to store
parms$store <- c(cond         = TRUE,
                 dat_full     = FALSE,
                 dat_miss     = FALSE,
                 sem_EST      = FALSE,
                 sem_CI       = FALSE,
                 lm_EST       = FALSE,
                 lm_CI        = FALSE,
                 miss_descrps = TRUE,
                 run_time_min = FALSE,
                 imp_values   = TRUE)

# Decide which methods you care about for the convergence check
parms$meth_sel <- data.frame(DURR_la    = TRUE,
                             DURR_SI    = FALSE,
                             IURR_la    = TRUE,
                             IURR_SI    = FALSE,
                             blasso     = TRUE,
                             blasso_SI  = FALSE,
                             bridge     = TRUE,
                             bridge_SI  = FALSE,
                             MI_PCA     = TRUE, # always true if SI is required
                             MI_CART    = TRUE,
                             MI_CART_SI = FALSE,
                             MI_RF      = TRUE,
                             MI_RF_SI   = FALSE,
                             MI_OP   = TRUE,
                             missFor = TRUE,
                             GS      = TRUE,
                             CC      = TRUE
)

parms$methods <- names(parms$meth_sel)[which(parms$meth_sel==TRUE)]

# Decide which paramters suit best the convergence check
parms$dt_rep     <- 2  # 10 # replications for averaging results (200 goal)
parms$chains     <- 2  # 5  # number of parallel chains for convergence check
parms$iters      <- 10 # 250 
parms$burnin_imp <- 0  # no need, I want to see the first part as well
parms$ndt        <- parms$chains # same chians as dataset in mice like imp
parms$thin       <- (parms$iters - parms$burnin_imp)/parms$ndt
parms$pos_dt     <- (parms$burnin_imp+1):parms$iters # candidate datasets (after convergence)
parms$keep_dt    <- parms$pos_dt[seq(1, length(parms$pos_dt), parms$thin)] # keep 1 dataset every thin

# For blasso
parms$chains_bl     <- 2  # 5 number of parallel chains for convergence check
parms$iters_bl      <- 10 # 2e3
parms$burnin_imp_bl <- 0  # how many imputation iterations should be discarded
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
conds <- conds[8, ]
  # condition 12 is the complex: interactions everywhere and larger dataset to "scan"
