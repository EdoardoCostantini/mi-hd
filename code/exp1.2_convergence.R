# Project:   imputeHD-comp
# Objective: Check convergence for collinearity study
# Author:    Edoardo Costantini
# Created:   2023-03-28
# Modified:  2023-03-28
# Notes: 

# 1. Convergence for all R-based MI methods -------------------------------------

# Clean environment
rm(list = ls())

# Source scripts
source("./init_general.R")
source("./exp1.2_init.R")

# 1.2 Modify parms and conds objects -------------------------------------------

# 1.2.1 Parameters -------------------------------------------------------------

# Which imputation method to use
parms$meth_sel <- data.frame(
    DURR_la = TRUE,
    IURR_la = TRUE,
    blasso = TRUE,
    bridge = TRUE,
    MI_PCA = TRUE,
    MI_CART = TRUE,
    MI_RF = TRUE,
    stepFor = FALSE,  # <- Special method requiring a different strategy
    MI_qp = TRUE,
    MI_am = TRUE,
    MI_OP = TRUE,
    missFor = TRUE,
    mean = TRUE,
    CC = TRUE,
    GS = TRUE
)
parms$methods <- names(parms$meth_sel)[which(parms$meth_sel == TRUE)]

# What to store
parms$store <- c(
    cond = TRUE,
    dat_full = TRUE,
    dat_miss = TRUE,
    sem_EST = FALSE,
    sem_CI = FALSE,
    lm_EST = FALSE,
    lm_CI = FALSE,
    fmi = FALSE,
    miss_descrps = FALSE,
    run_time_min = FALSE,
    imp_values = TRUE # <- this is what we want for traceplots
)

# Iterations, repetitions, etc
parms$dt_rep <- 10
parms$ndt <- parms$mice_ndt <- parms$chains <- 5
parms$iters <- parms$mice_iters <- 10 # 250 iterations
parms$burnin_imp <- 0
parms$thin <- (parms$iters - parms$burnin_imp) / parms$ndt
parms$pos_dt <- (parms$burnin_imp + 1):parms$iters # candidate datasets (after convergence)
parms$keep_dt <- parms$pos_dt[seq(
    1,
    length(parms$pos_dt),
    parms$thin
)] # keep 1 dataset every thin

# Report names
filename <- paste0(
    "exp",
    parms$exp,
    "_convergence_all_meth_",
    parms$start_time
)
parms$report_file_name <- paste0(filename, ".txt")
parms$results_file_name <- paste0(filename, ".rds")

# 1.2.2 Conditions -------------------------------------------------------------

latent <- FALSE
pm <- .3
p <- 500 # c(50, 500) # number of variables
collinearity <- c(.99)
ridge <- 0.01

# Create experimental conditions
conds <- expand.grid(
    ridge = ridge,
    p = p,
    latent = latent,
    pm = pm,
    collinearity = collinearity
)

# Bridge special parameters per condition
conds$ridge <- 0.01

# 1.3 Perform simulation -------------------------------------------------------

# Start simulation 
sim_start <- Sys.time()

# Create a cluster object:
clus <- makeCluster(parms$dt_rep)

# Export call for packages
clusterEvalQ(cl = clus, expr = source("./init_general.R"))

# Export parms and conds objects
clusterExport(
    cl = clus,
    varlist = list(
        "parms",
        "conds"
    )
)

# Run the computations in parallel on the 'clus' object:
out <- parLapply(
    cl = clus,
    X = 1:(parms$dt_rep),
    fun = doRep,
    conds = conds,
    parms = parms
)

# Kill the cluster:
stopCluster(clus)

sim_ends <- Sys.time()

# Append Parms and conds to out object
out$parms <- parms
out$conds <- conds
out$session_info <- devtools::session_info()

# Save results
saveRDS(
    out,
    paste0(
        parms$outDir,
        parms$results_file_name
    )
)

# 1.4 Analyze results ----------------------------------------------------------

# Load data

# 20200803 Correct
out_cnv <- readRDS("")

# What to show
iters_range <- 1:5 # which set of iterations
y_range <- c(3, 6)
exp_dat <- 1 # which data replication (10 possibilities)

# Run description
data.frame(
    Chians = c(all = out_cnv$parms$chains, MICE = out_cnv$parms$mice_ndt),
    Iters = c(all = out_cnv$parms$iters, MICE = out_cnv$parms$mice_iters)
)

# DURR_la
mean_traceplot(
    out = out_cnv,
    method = out_cnv$parms$method[1],
    dat = exp_dat,
    y_range = y_range,
    iters = iters_range
)

# IURR_la
mean_traceplot(out_cnv,
    method = out_cnv$parms$method[2],
    dat = exp_dat,
    y_range = y_range, iters = iters_range
)
# bridge
mean_traceplot(out_cnv,
    method = out_cnv$parms$method[3],
    dat = exp_dat,
    y_range = c(3, 7), iters = iters_range
)
# blasso
mean_traceplot(out_cnv,
    method = out_cnv$parms$method[4],
    dat = exp_dat,
    y_range = y_range, iters = iters_range
)
# MI_PCA
mean_traceplot(out_cnv,
    method = out_cnv$parms$method[5],
    dat = exp_dat,
    y_range = y_range, iters = iters_range
)
# MI_CART
mean_traceplot(out_cnv,
    method = out_cnv$parms$method[6],
    dat = exp_dat,
    y_range = y_range, iters = iters_range
)
# MI_RF
mean_traceplot(out_cnv,
    method = out_cnv$parms$method[7],
    dat = exp_dat,
    y_range = y_range, iters = iters_range
)
# MI_T
mean_traceplot(out_cnv,
    method = out_cnv$parms$method[8],
    dat = exp_dat,
    y_range = y_range, iters = iters_range
)
# MI_am
mean_traceplot(out_cnv,
    method = out_cnv$parms$method[8],
    dat = exp_dat,
    y_range = y_range, iters = iters_range
)
# MI_OP
mean_traceplot(out_cnv,
    method = out_cnv$parms$method[2],
    dat = exp_dat,
    y_range = y_range,
    iters = iters_range
)

# 2. Convergence for IVEware ---------------------------------------------------

# Clean environment
rm(list = ls())

# Source scripts
source("./init_general.R")
source("./exp1.2_init.R")

# 2.2 Modify parms and conds objects -------------------------------------------

# 2.2.1 Parameters -------------------------------------------------------------

# Which imputation method to use
parms$meth_sel <- data.frame(
    DURR_la = FALSE,
    IURR_la = FALSE,
    blasso = FALSE,
    bridge = FALSE,
    MI_PCA = FALSE,
    MI_CART = FALSE,
    MI_RF = FALSE,
    stepFor = TRUE, # <- method of interest
    MI_qp = FALSE,
    MI_am = FALSE,
    MI_OP = FALSE,
    missFor = TRUE,
    mean = TRUE,
    CC = TRUE,
    GS = TRUE
)
parms$methods <- names(parms$meth_sel)[which(parms$meth_sel == TRUE)]

# What to store
parms$store <- c(
    cond = TRUE,
    dat_full = TRUE,
    dat_miss = TRUE,
    sem_EST = FALSE,
    sem_CI = FALSE,
    lm_EST = FALSE,
    lm_CI = FALSE,
    fmi = FALSE, #  <- this is what we want for CV
    miss_descrps = FALSE,
    run_time_min = FALSE,
    imp_values = TRUE
)

# Iterations, repetitions, etc
parms$dt_rep <- 3
parms$ndt <- 5
parms$iters <- 5
parms$burnin_imp <- 1
parms$thin <- (parms$iters - parms$burnin_imp) / parms$ndt
parms$pos_dt <- (parms$burnin_imp + 1):parms$iters # candidate datasets (after convergence)
parms$keep_dt <- parms$pos_dt[seq(
    1,
    length(parms$pos_dt),
    parms$thin
)] # keep 1 dataset every thin

# Report names
filename <- paste0(
    "exp",
    parms$exp, "_",
    "cv_IVEware_",
    parms$start_time
)
parms$report_file_name <- paste0(filename, ".txt")
parms$results_file_name <- paste0(filename, ".rds")

# 2.2.2 Conds ------------------------------------------------------------------

latent <- FALSE
pm <- .3
p <- c(50, 500) # c(50, 500) # number of variables
collinearity <- c(.6, .8, .90, .99)
minR2 <- c(10^seq(from = -1, to = -7, by = -1))

# Create experimental conditions
conds <- expand.grid(
    minR2 = minR2,
    p = p,
    latent = latent,
    pm = pm,
    collinearity = collinearity
)[1:4, ]

# IVEware special parameters per condition #TODO: cross-validate
conds$ridge <- NA

# 2.3 Perform simulation -------------------------------------------------------

# Start simulation
sim_start <- Sys.time()

# Create a cluster object:
clus <- makeCluster(parms$dt_rep)

# Export call for packages
clusterEvalQ(cl = clus, expr = source("./init_general.R"))

# Export parms and conds objects
clusterExport(
    cl = clus,
    varlist = list(
        "parms",
        "conds"
    )
)

# Run the computations in parallel on the 'clus' object:
out <- parLapply(
    cl = clus,
    X = 1:(parms$dt_rep),
    fun = doRep,
    conds = conds,
    parms = parms
)

# Kill the cluster:
stopCluster(clus)

sim_ends <- Sys.time()

# Append Parms and conds to out object
out$parms <- parms
out$conds <- conds
out$session_info <- devtools::session_info()

# Save results
saveRDS(
    out,
    paste0(
        parms$outDir,
        parms$results_file_name
    )
)

# 2.4 Analyze results ----------------------------------------------------------

# Load data
out <- readRDS("../output/exp1_cv_bridge_20220224_1042.rds")

# Obtain conditions with cv ridge
conds_bridge <- cvParm(
    out,
    cv.parm = "ridge",
    exp_factors = colnames(out$conds)[c(2, 4)]
)

# Selected values
conds_bridge$values

# Plot
conds_bridge$plot
