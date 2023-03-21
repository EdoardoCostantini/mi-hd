# Title:    Simulation script for ridge paramter cross-validation
# Project:  Imputing High Dimensional Data
# Author:   Edoardo Costantini
# Created:  2020-08-25
# Modified: 2022-02-23

rm(list = ls())
source("./init_general.R")
source("./exp1_init.R")

# Cross-validation -------------------------------------------------------------

# > Modify parameter of interest -----------------------------------------------

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

# What to store (FMI for cross-validation)
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

# > Conditions -----------------------------------------------------------------

minR2 <- c(.1, .01)
p <- c(50, 500) # c(50, 500) # number of variables
latent <- c(FALSE, TRUE)[1]
pm <- c(.1, .3)

conds <- expand.grid(
        minR2 = minR2,
        p = p,
        latent = latent,
        pm = pm
)

# > Run ------------------------------------------------------------------------

# Define empty list
out <- vector("list", parms$dt_rep)

# Loop over repetitions
for (r in 1:parms$dt_rep) {
        out[[r]] <- doRep(r, conds = conds, parms = parms)
}

# Attach extra objects
out$parms <- parms
out$conds <- conds
out$session_info <- devtools::session_info()

# Save object
saveRDS(
        out,
        paste0(
                parms$outDir,
                parms$results_file_name
        )
)

# Obtain conditions with cv ridge
conds_bridge <- cvParm(
        out = out,
        cv.parm = colnames(conds)[1],
        exp_factors = colnames(conds)[c(2, 4)]
)

# Look at plot
conds_bridge$plot

# Look at solutions
conds_bridge$values

# Running simulation -----------------------------------------------------------
