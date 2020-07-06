### Title:    Convergence Checks: Initialization Script for Exp 1 
### Project:  Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-07-03
### Notes:    The goal is to find out how many iterations should be used
###           in the full study for each method.

# Get uninteresting parameters from the init_exp1.R script
source("./init_exp1.R")

# Decide which methods you car about for the convergence check
parms$meth_sel <- data.frame(DURR_la = TRUE,
                             DURR_el = FALSE,
                             IURR_la = TRUE,
                             IURR_el = FALSE,
                             blasso  = TRUE,
                             bridge  = TRUE,
                             MI_PCA  = TRUE,
                             MI_CART = TRUE,
                             MI_RF   = TRUE,
                             MI_OP   = TRUE,
                             missFor = TRUE,
                             GS      = TRUE,
                             CC      = TRUE
)

# Decide which paramters suit best the convergence check
parms$dt_rep     <- 10 # replications for averaging results (200 goal)
parms$chains     <- 5  # number of parallel chains for convergence check
parms$iters      <- 200 
parms$burnin_imp <- 0  # no need, I want to see the first part as well
parms$ndt        <- parms$chains # same chians as dataset in mice like imp
parms$thin       <- (parms$iters - parms$burnin_imp)/parms$ndt
parms$pos_dt     <- (parms$burnin_imp+1):parms$iters # candidate datasets (after convergence)
parms$keep_dt    <- parms$pos_dt[seq(1, length(parms$pos_dt), parms$thin)] # keep 1 dataset every thin

# For mice-like algorithms
parms$mice_iters <- parms$iters # for this they are the same.
parms$mice_ndt   <- parms$ndt # 10 # number of imputed datasets to pool esitmaes from (10)

# Change name of outputs

# Output and Progres report related
parms$outDir <- "../output/"
parms$start_time <- format(Sys.time(), "%Y%m%d_%H%M")
parms$report_file_name <- paste0("cnv_check_", 
                                 parms$start_time, 
                                 ".txt")
parms$results_file_name <- paste0("cnv_check_", 
                                  parms$start_time,
                                  ".rds")
parms$description <- c("For the most challanging condition, data generation is repeated 
                        a few times, imputation models are used. Imputations are can 
                        then be used to check convergence")

# Fix condition to desired one for convergence check
p   <- 500
pm <- .3
conds <- data.frame(p = p, pm = pm)

