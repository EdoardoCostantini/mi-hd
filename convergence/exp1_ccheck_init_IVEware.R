# Title:    Convergence Checks: Initialization Script for IVEware Exp 1
# Project:  Imputing High Dimensional Data
# Author:   Edoardo Costantini
# Created:  2023-03-24
# Modified: 2023-03-27
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
    stepFor = TRUE, # The only MI method that should be active for this script
    MI_qp = FALSE,
    MI_am = FALSE,
    MI_OP = FALSE,
    missFor = TRUE,
    mean = TRUE,
    CC = TRUE,
    GS = TRUE
)

# Define the number of multiply imputed datasets
parms$ndt <- 30

# Output and progress report related
parms$outDir <- "../output/"
parms$start_time <- format(Sys.time(), "%Y%m%d_%H%M")
parms$report_file_name <- paste0(
    "exp", parms$exp,
    "_cv_IVEware_",
    parms$start_time,
    ".txt"
)
parms$results_file_name <- paste0(
    "exp", parms$exp,
    "_cv_IVEware_",
    parms$start_time,
    ".rds"
)
parms$description <- c(
    "For the most challenging condition, one data is generated and, IVEware imputation is used. The final imputations can be compared with the observed densities."
)

# Fix condition to desired one for convergence check
p <- 500
pm <- .3
iters <- c(1, 5, 10, 20, 40, 80, 160, 240, 320)
seed <- c(2023, 3202, 1610, 4004, 9876)
conds <- expand.grid(p = p, pm = pm, iters = iters, seed = seed)
