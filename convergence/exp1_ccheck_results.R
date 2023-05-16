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
m <- 5
# Number of 

# Create an empty object for storing the results of the repetition
mids_long_r <- NULL

for(r in 1:reps){

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
    mids_long_m <- NULL

    # For every method
    for (m in 1:length(out_cnv[[r]][[1]]$imp_values)) {

        # Create an empty object for storing the results of chain
        mids_long_cn <- NULL

        # Is the object mids?
        is.mids <- class(out_cnv[[r]][[1]]$imp_values[[m]]) == "mids"

        # If so, do this
        if(is.mids == FALSE){
            # For every chain
            for (cn in 1:length(out_cnv[[r]][[1]]$imp_values[[m]])) {

                # Create an empty object for storing the results for the variable
                mids_long_v <- NULL

                # For every variable
                for (v in 1:length(out_cnv[[r]][[1]]$imp_values[[m]][[cn]])) {
                        # Create a unit of information
                        temp_mean_sd <- data.frame(
                            repetition = r,
                            method = names(out_cnv[[r]][[1]]$imp_values)[m],
                            chain = cn,
                            iter = 1:nrow(out_cnv[[r]][[1]]$imp_values[[m]][[cn]][[v]]),
                            variable = paste0("z", c(1:3, 6:8)[v]),
                            mean = rowMeans(out_cnv[[r]][[1]]$imp_values[[m]][[cn]][[v]]),
                            sd = apply(out_cnv[[r]][[1]]$imp_values[[m]][[cn]][[v]], 1, sd)
                        )
                    # Reshape for future plots
                    mids_long <- reshape2::melt(
                        temp_mean_sd,
                        id.var = c("repetition", "method", "chain", "iter", "variable"),
                        variable.name = "outcome"
                    )

                    # Combine different variables
                    mids_long_v <- rbind(mids_long_v, mids_long)
                }
                # Combine different chains
                mids_long_cn <- rbind(mids_long_cn, mids_long_v)
            }
        } else {
            # Obtain an object on the same shame starting from
            id_values <- rowSums(is.nan(out_cnv[[r]][[1]]$imp_values[[m]]$chainMean)) == 0

            # Keep only values
            mids_chainMean <- out_cnv[[r]][[1]]$imp_values[[m]]$chainMean[id_values, ,]
            mids_chainVar <- out_cnv[[r]][[1]]$imp_values[[m]]$chainVar[id_values, ,]
        }

        # Combine different methods
        mids_long_m <- rbind(mids_long_m, mids_long_cn)
    }

    # Combine repetitions
    mids_long_r <- rbind(mids_long_r, mids_long_m)
}