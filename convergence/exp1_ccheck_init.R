# Project:   imputeHD-comp
# Objective: Convergence Checks: Initialization Script for Exp 1
# Author:    Edoardo Costantini
# Created:   2020-07-03
# Modified:  2023-03-28
# Notes:

# Get uninteresting parameters from the init_exp1.R script
source("./exp1_init.R")

# Decide which methods you car about for the convergence check
parms$meth_sel <- data.frame(
    DURR_la = TRUE,
    IURR_la = TRUE,
    blasso = TRUE,
    bridge = TRUE,
    MI_PCA = TRUE,
    MI_CART = TRUE,
    MI_RF = TRUE,
    stepFor = FALSE,
    MI_qp = TRUE,
    MI_am = TRUE,
    MI_OP = TRUE,
    missFor = TRUE,
    mean = TRUE,
    CC = TRUE,
    GS = TRUE
)

# Decide which parameters suit best the convergence check
parms$dt_rep <- 10 # replications for averaging results (200 goal)
parms$chains <- 5 # number of parallel chains for convergence check
parms$iters <- 250
parms$burnin_imp <- 0 # no need, I want to see the first part as well
parms$ndt <- parms$chains # same chains as dataset in mice like imp
parms$thin <- (parms$iters - parms$burnin_imp) / parms$ndt
parms$pos_dt <- (parms$burnin_imp + 1):parms$iters # candidate datasets (after convergence)
parms$keep_dt <- parms$pos_dt[seq(1, length(parms$pos_dt), parms$thin)] # keep 1 dataset every thin

# For mice-like algorithms
parms$mice_iters <- parms$iters # for this they are the same.
parms$mice_ndt <- parms$ndt # 10 # number of imputed datasets to pool esitmaes from (10)

# Change name of outputs

# Output and Progres report related
parms$outDir <- "../output/"
parms$start_time <- format(Sys.time(), "%Y%m%d_%H%M")
parms$report_file_name <- paste0(
    "exp1_cnv_check_",
    parms$start_time,
    ".txt"
)
parms$results_file_name <- paste0(
    "exp1_cnv_check_",
    parms$start_time,
    ".rds"
)
parms$description <- c("For the most challanging condition, data generation is repeated
                        a few times, imputation models are used. Imputations are can
                        then be used to check convergence")

# Fix condition to desired one for convergence check
latent <- FALSE
p <- 500
pm <- .3
collinearity <- c(.99)
ridge <- .01
minR2 <- 1e-3

# Create experimental conditions
conds <- expand.grid(
    p = p,
    latent = latent,
    pm = pm,
    collinearity = collinearity
)

# Bridge special parameters per condition
conds$ridge <- ridge

# IVEware special parameters per condition
conds$minR2 <- minR2
