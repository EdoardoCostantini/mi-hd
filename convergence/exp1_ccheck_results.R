# Project:   imputeHD-comp
# Objective: Convergence checks for all methods in EXP 1 Evs
# Author:    Edoardo Costantini
# Created:   2020-07-03
# Modified:  2023-03-28
# Notes:

rm(list = ls())
source("./init_general.R")

# Read Simulation Results -------------------------------------------------
# Special Run with multiple chains. This output does not contain results
# of the simualtion study. It contains an imputation run with many iterations,
# many chians, for the most difficult condition of experiment 1.

# 20200803 Correct
out_cnv <- readRDS("../output/exp1_conv_20200731_1652.rds")

# What to show
iters_range <- 1:250 # which set of iterations
y_range <- c(3, 6)
exp_dat <- 1 # which data replication (10 possibilities)

# Run description
data.frame(
    Chians = c(all = out_cnv$parms$chains, MICE = out_cnv$parms$mice_ndt),
    Iters = c(all = out_cnv$parms$iters, MICE = out_cnv$parms$mice_iters)
)

# DURR_la
mean_traceplot(out_cnv,
    method = out_cnv$parms$method[1],
    dat = exp_dat,
    y_range = y_range, iters = iters_range
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
mean_traceplot(
    out = out_cnv,
    method = out_cnv$parms$method[5],
    dat = 1,
    y_range = y_range,
    iters = iters_range
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
mean_traceplot(
    out = out_cnv,
    method = out_cnv$parms$method[8],
    dat = exp_dat,
    y_range = y_range,
    iters = iters_range
)

# Resize object ----------------------------------------------------------------

# Number of repetitions
reps <- length(out_cnv) - 1
r <- 1
m <- 1
v <- 1

# Create an empty object for storing the results of the repetition
mids_r <- NULL

for (r in 1:reps) {
    # Drop unused methods
    out_cnv[[r]][[1]]$imp_values <- out_cnv[[r]][[1]]$imp_values[!grepl("_el", names(out_cnv[[r]][[1]]$imp_values))]

    # Rename methods
    names(out_cnv[[r]][[1]]$imp_values) <- gsub("_la", "", names(out_cnv[[r]][[1]]$imp_values))
    names(out_cnv[[r]][[1]]$imp_values) <- gsub("bridge", "BRidge", names(out_cnv[[r]][[1]]$imp_values))
    names(out_cnv[[r]][[1]]$imp_values) <- gsub("blasso", "BLasso", names(out_cnv[[r]][[1]]$imp_values))
    names(out_cnv[[r]][[1]]$imp_values) <- gsub("_", "-", names(out_cnv[[r]][[1]]$imp_values))

    # Create a storing object
    mids_methods <- vector("list", length(out_cnv[[r]][[1]]$imp_values))
    names(mids_methods) <- names(out_cnv[[r]][[1]]$imp_values)

    # Create an empty object for storing the results of the method
    mids_method <- NULL

    # For every method
    for (m in 1:length(out_cnv[[r]][[1]]$imp_values)) {
    
        # Is the object mids?
        is.mids <- class(out_cnv[[r]][[1]]$imp_values[[m]]) == "mids"

        # If so, do this
        if (is.mids == FALSE) {
            # Create two empty arrays
            mids_chainMean <- array(
                NA,
                dim = c(
                    length(out_cnv$parms$z_m_id), # number imputed vars
                    out_cnv$parms$iters, # number of iterations
                    out_cnv$parms$ndt # number of chains
                ),
                dimnames = list(
                    paste0("z", out_cnv$parms$z_m_id),
                    1:out_cnv$parms$iters,
                    paste0("chain", 1:out_cnv$parms$ndt)
                )
            )
            mids_chainVar <- array(
                NA,
                dim = c(
                    length(out_cnv$parms$z_m_id), # number imputed vars
                    out_cnv$parms$iters, # number of iterations
                    out_cnv$parms$ndt # number of chains
                ),
                dimnames = list(
                    paste0("z", out_cnv$parms$z_m_id),
                    1:out_cnv$parms$iters,
                    paste0("chain", 1:out_cnv$parms$ndt)
                )
            )

            # For every chain
            for (cn in 1:length(out_cnv[[r]][[1]]$imp_values[[m]])) {
                # For every variable
                for (v in 1:length(out_cnv[[r]][[1]]$imp_values[[m]][[cn]])) {
                    # Compute mean imputed values
                    mids_chainMean[v, , cn] <- rowMeans(out_cnv[[r]][[1]]$imp_values[[m]][[cn]][[v]])

                    # Compute sd of imputed values
                    mids_chainVar[v, , cn] <- apply(out_cnv[[r]][[1]]$imp_values[[m]][[cn]][[v]], 1, var)
                }

                # Create slim mids objects
                mids_m <- list(
                    chainMean = mids_chainMean,
                    chainVar = mids_chainVar,
                    m = length(out_cnv[[r]][[1]]$imp_values[[m]]),
                    iteration = ncol(mids_chainMean)
                )
            }
        } else {
            # Obtain an object on the same shame starting from
            id_values <- rowSums(is.nan(out_cnv[[r]][[1]]$imp_values[[m]]$chainMean)) == 0

            # Create slim mids objects
            mids_m <- list(
                chainMean = out_cnv[[r]][[1]]$imp_values[[m]]$chainMean[id_values, , ],
                chainVar = out_cnv[[r]][[1]]$imp_values[[m]]$chainVar[id_values, , ],
                m = out_cnv[[r]][[1]]$imp_values[[m]]$m,
                iteration = out_cnv[[r]][[1]]$imp_values[[m]]$iteration
            )
        }

        mids_method[[m]] <- mids_m

    }

    # Give names to the mids objects
    names(mids_method) <- names(out_cnv[[r]][[1]]$imp_values)

    # Combine repetitions
    mids_r[[r]] <- mids_method
}

# Save slim and prepped mids object
out_cnv <- saveRDS(mids_r, "../output/exp1_conv_20200731_1652_shiny.rds")