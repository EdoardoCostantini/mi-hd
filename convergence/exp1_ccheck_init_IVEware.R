# Title:    Convergence Checks: Initialization Script for IVEware Exp 1
# Project:  Imputing High Dimensional Data
# Author:   Edoardo Costantini
# Created:  2023-03-24
# Modified: 2023-03-24
# Notes:    The goal is to find out how many iterations should be used
#           in the full study for each method.

# Decide which methods you car about for the convergence check
parms$meth_sel <- data.frame(
    DURR_la = FALSE,
    IURR_la = FALSE,
    blasso = FALSE,
    bridge = FALSE,
    MI_PCA = FALSE,
    MI_CART = FALSE,
    MI_RF = FALSE,
    stepFor = TRUE,
    MI_qp = FALSE,
    MI_am = FALSE,
    MI_OP = FALSE,
    missFor = TRUE,
    mean = TRUE,
    CC = TRUE,
    GS = TRUE
)

# Decide which parameters suit best the convergence check
parms$dt_rep <- 1 # replications for averaging results (200 goal)
parms$chains <- 10 # number of parallel chains for convergence check
parms$iters <- 250
parms$burnin_imp <- 0 # no need, I want to see the first part as well
parms$ndt <- parms$chains # same chians as dataset in mice like imp
parms$thin <- (parms$iters - parms$burnin_imp) / parms$ndt
parms$pos_dt <- (parms$burnin_imp + 1):parms$iters # candidate datasets (after convergence)
parms$keep_dt <- parms$pos_dt[seq(1, length(parms$pos_dt), parms$thin)] # keep 1 dataset every thin

# For mice-like algorithms
parms$mice_iters <- parms$iters # for this they are the same.
parms$mice_ndt <- parms$ndt # 10 # number of imputed datasets to pool esitmaes from (10)

# Output and progress report related
parms$outDir <- "../output/"
parms$start_time <- format(Sys.time(), "%Y%m%d_%H%M")
parms$report_file_name <- paste0(
    "cnv_check_",
    parms$start_time,
    ".txt"
)
parms$results_file_name <- paste0(
    "cnv_check_",
    parms$start_time,
    ".rds"
)
parms$description <- c(
    "For the most challenging condition, one data is generated and, IVEware imputation is used. The final imputations can be compared with the observed densities."
)

# Fix condition to desired one for convergence check
p <- 500
pm <- .3
iters <- c(1, 5, 10, 15, 20, 30, 70, 150, 250)
conds <- data.frame(p = p, pm = pm, iters = iters)
