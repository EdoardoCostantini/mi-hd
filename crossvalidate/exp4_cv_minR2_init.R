# Project:   imputeHD-comp
# Objective: Edits to default parms for cv of minR2 in IVEware imputation for EVS
# Author:    Edoardo Costantini
# Created:   2023-03-22
# Modified:  2023-03-22
# Notes:

# Modify parameter of interest -------------------------------------------------

# Which imputation method to use
parms$meth_sel <- data.frame(
    DURR_la = FALSE,
    IURR_la = FALSE,
    blasso = FALSE,
    bridge = TRUE,
    MI_PCA = FALSE,
    MI_CART = FALSE,
    MI_RF = FALSE,
    MI_qp = FALSE,
    MI_am = FALSE,
    MI_OP = FALSE,
    missFor = TRUE,
    mean = TRUE,
    CC = TRUE,
    GS = TRUE
)

parms$methods <- names(parms$meth_sel)[which(parms$meth_sel == TRUE)]

# Store only FMI for cross-validation
parms$store <- c(
    cond = TRUE,
    dat_full = FALSE,
    dat_miss = FALSE,
    sem_EST = FALSE,
    sem_CI = FALSE,
    lm_EST = FALSE,
    lm_CI = FALSE,
    fmi = TRUE,
    miss_descrps = FALSE,
    run_time_min = FALSE,
    imp_values = FALSE
)

# Iterations, repetitions, etc
parms$dt_rep <- 20
parms$chains <- 1
parms$iters <- 70
parms$burnin_imp <- 50
parms$ndt <- 10
parms$thin <- (parms$iters - parms$burnin_imp) / parms$ndt
parms$pos_dt <- (parms$burnin_imp + 1):parms$iters
parms$keep_dt <- parms$pos_dt[seq(1, length(parms$pos_dt), parms$thin)]

# Report names
parms$report_file_name <- paste0(
    "exp",
    parms$exp, "_",
    "cv_IVEware_",
    parms$start_time,
    ".txt"
)
parms$results_file_name <- paste0(
    "exp",
    parms$exp, "_",
    "cv_IVEware_",
    parms$start_time,
    ".rds"
)

# Conditions -------------------------------------------------------------------

# Experimental factors
n <- c(1e3, 3e2) # number of observations

# Dataframe of conditions w/ ridge
minR2 <- c(.2, 10^seq(from = -1, to = -5, by = -1))

# Define experimental conditions
conds <- expand.grid(
    minR2 = minR2,
    n = n
)
