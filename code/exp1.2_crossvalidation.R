# Project:   imputeHD-comp
# Objective: Script to cross-validate ridge and minR2 for collinearity study
# Author:    Edoardo Costantini
# Created:   2023-03-28
# Modified:  2023-04-08
# Notes: 

# 1. Cross-validation of ridge in Bridge ---------------------------------------

# Clean environment
rm(list = ls())

# Source scripts
source("./init_general.R")
source("./exp1.2_init.R")

# 1.2 Modify parms and conds objects -------------------------------------------

# 1.2.1 Parameters -------------------------------------------------------------

# Which imputation method to use
parms$meth_sel <- data.frame(
    DURR_la = FALSE,
    IURR_la = FALSE,
    blasso = FALSE,
    bridge = TRUE, # <- method of interest
    MI_PCA = FALSE,
    MI_CART = FALSE,
    MI_RF = FALSE,
    stepFor = FALSE,
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
    dat_full = FALSE,
    dat_miss = FALSE,
    sem_EST = FALSE,
    sem_CI = FALSE,
    lm_EST = FALSE,
    lm_CI = FALSE,
    fmi = TRUE,            #  <- this is what we want for CV
    miss_descrps = FALSE,
    run_time_min = FALSE,
    imp_values = FALSE
)

# Iterations, repetitions, etc
parms$dt_rep <- 20
parms$ndt <- 10
parms$iters <- 70
parms$burnin_imp <- 50
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
    "cv_bridge_",
    parms$start_time
)
parms$report_file_name <- paste0(filename, ".txt")
parms$results_file_name <- paste0(filename, ".rds")

# 1.2.2 Conditions -------------------------------------------------------------

latent <- FALSE
pm <- .3
p <- c(50, 500) # c(50, 500) # number of variables
collinearity <- c(.6, .8, .90)
ridge <- 10^seq(from = -1, to = -8, by = -1)

# Create experimental conditions
conds <- expand.grid(
    ridge = ridge,
    p = p,
    latent = latent,
    pm = pm,
    collinearity = collinearity
)

# IVEware special parameters per condition #TODO: cross-validate
conds$minR2 <- NA

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
out <- readRDS("../output/exp1_2_cv_bridge_20230405_1449.rds")

# Obtain conditions with cv ridge
conds_bridge <- cvParm(
    out,
    cv.parm = "ridge",
    mods = "sems",
    exp_factors = colnames(out$conds)[c(2, 5)]
)

# Selected values
conds_bridge$solution
conds_bridge$solution_1se

# Plot
conds_bridge$plot

# 2. Cross-validation of minR2 in IVEware --------------------------------------

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
    dat_full = FALSE,
    dat_miss = FALSE,
    sem_EST = FALSE,
    sem_CI = FALSE,
    lm_EST = FALSE,
    lm_CI = FALSE,
    fmi = TRUE, #  <- this is what we want for CV
    miss_descrps = FALSE,
    run_time_min = FALSE,
    imp_values = FALSE
)

# Iterations, repetitions, etc
parms$dt_rep <- 30
parms$ndt <- 50
parms$iters <- 50
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
collinearity <- c(.6, .8, .90)
minR2 <- c(10^seq(from = -1, to = -7, by = -1))

# Create experimental conditions
conds <- expand.grid(
    minR2 = minR2,
    p = p,
    latent = latent,
    pm = pm,
    collinearity = collinearity
)

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
out <- readRDS("../output/exp1_2_cv_IVEware_20230406_1053.rds")

# Obtain conditions with cv ridge
conds_minR2 <- cvParm(
    out,
    mods = names(out[[1]][[1]]$fmi)[1],
    cv.parm = "minR2",
    exp_factors = colnames(out$conds)[c(2, 5)]
)

# Selected values
conds_minR2$solution
conds_minR2$solution_1se

# Plot
conds_minR2$plot
