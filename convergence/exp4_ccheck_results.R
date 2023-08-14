# Project:   imputeHD-comp
# Objective: Convergence checks for all methods in EXP 4 Evs
# Author:    Edoardo Costantini
# Created:   2020-07-03
# Modified:  2023-03-28
# Notes: 

rm(list = ls())
source("./init_general.R")

# Read Simulation Results -------------------------------------------------
# Special Run with multiple chains. This output does not contain results
# of the simualtion study. It contains an imputation run with many iterations,
# many chians, for the most difficult condition of experiment 3.

out_cnv <- readRDS("../output/exp4_conv_20201015_1555.rds")

# Run description

data.frame(
  Reps   = c(out_cnv$parms$dt_rep, out_cnv$parms$dt_rep, out_cnv$parms$dt_rep),
  Chians = c(all = out_cnv$parms$chains, 
             blasso = out_cnv$parms$chains,
             MICE = out_cnv$parms$mice_ndt),
  Iters = c(all = out_cnv$parms$iters, 
            blasso = out_cnv$parms$iters_bl,
            MICE = out_cnv$parms$mice_iters)
)

# R-hat -------------------------------------------------------------------

# Perform for all methods that support it

  # Number of dataset for which I have convergence checks
  out_cnv$parms$dt_rep

  # Methods for which I have R-hat
  mth_range <- c(1:4, 6:7)
  names(mth_range) <- out_cnv$parms$method[mth_range]

  # Average Rhats (over data repetitions)
  Rhats <- sapply(mth_range, function(m){
    rowMeans(sapply(1:out_cnv$parms$dt_rep, 
                    Rhat.sim,
                    out  = out_cnv,
                    cond = 1,
                    meth = out_cnv$parms$method[m],
                    iter_max = out_cnv$parms$iters,
                    iter_burn = 50))
  })

# Plots -------------------------------------------------------------------

  # What to show
  exp_dat     <- 2
  iters_range <- 1:250 # which set of iterations
  y_range     <- c(2, 3.5)

  # Print ALL
  lapply(seq_along(out_cnv$parms$method[1:8]), function(x){
    mean_traceplot(out_cnv, 
                   method = out_cnv$parms$method[x], 
                   dat = exp_dat, 
                   y_range = y_range, 
                   iters = iters_range)}
  )

# Study specific imputations
# DURR_la (50-100 good for all)
print(out_cnv$parms$method[1])
mean_traceplot(out_cnv, 
               method = out_cnv$parms$method[1], 
               dat = exp_dat, 
               y_range = y_range, 
               iters = iters_range)

# IURR_la (50-100 good for all)
print(out_cnv$parms$method[2])
mean_traceplot(out_cnv, 
               method = out_cnv$parms$method[2], 
               dat = exp_dat, 
               y_range = y_range, 
               iters = iters_range)

# blasso (200-300 is a good range to chose from)
print(out_cnv$parms$method[3])
iters_range_bl <- 1:250
mean_traceplot(out_cnv,
               method   = out_cnv$parms$method[3],
               dat      = exp_dat, 
               y_range  = c(2, 3.5), iters = iters_range)

# bridge (50-100 good for all)
print(out_cnv$parms$method[4])
mean_traceplot(out_cnv, 
               method   = out_cnv$parms$method[4], 
               dat      = exp_dat,
               y_range  = c(-100, 100),
               iters    = iters_range)

# MI_PCA (50-100 good for all)
print(out_cnv$parms$method[5])
mean_traceplot(out_cnv,
               method = out_cnv$parms$method[5], 
               dat = exp_dat, 
               y_range = y_range, iters = iters_range)

# MI_CART
print(out_cnv$parms$method[6])
mean_traceplot(out_cnv,
               method = out_cnv$parms$method[6], 
               dat = exp_dat, 
               y_range = y_range, iters = iters_range)

# MI_RF
print(out_cnv$parms$method[7])
mean_traceplot(out_cnv,
               method = out_cnv$parms$method[7],
               dat = exp_dat, 
               y_range = y_range, iters = iters_range)

# MI_T
print(out_cnv$parms$method[8])
mean_traceplot(out_cnv, 
               method = out_cnv$parms$method[8], 
               dat = exp_dat, 
               y_range = y_range, iters = iters_range)

# 1.5 Reshape for shiny app ----------------------------------------------------

# Number of repetitions
reps <- length(out_cnv) - 2
r <- 1
m <- 1
v <- 1

# Create an empty object for storing the results of the repetition
mids_r <- NULL

for (r in 1:reps) {
  # Rename methods
  names(out_cnv[[r]][[1]]$imp_values) <- gsub("_la", "", names(out_cnv[[r]][[1]]$imp_values))
  names(out_cnv[[r]][[1]]$imp_values) <- gsub("bridge", "BRidge", names(out_cnv[[r]][[1]]$imp_values))
  names(out_cnv[[r]][[1]]$imp_values) <- gsub("blasso", "BLasso", names(out_cnv[[r]][[1]]$imp_values))
  names(out_cnv[[r]][[1]]$imp_values) <- gsub("_", "-", names(out_cnv[[r]][[1]]$imp_values))
  names(out_cnv[[r]][[1]]$imp_values) <- gsub("qp", "QP", names(out_cnv[[r]][[1]]$imp_values))
  names(out_cnv[[r]][[1]]$imp_values) <- gsub("am", "AM", names(out_cnv[[r]][[1]]$imp_values))
  names(out_cnv[[r]][[1]]$imp_values) <- gsub("OP", "OR", names(out_cnv[[r]][[1]]$imp_values))

  # Create a storing object
  mids_methods <- vector("list", length(out_cnv[[r]][[1]]$imp_values))
  names(mids_methods) <- names(out_cnv[[r]][[1]]$imp_values)

  # Create an empty object for storing the results of the method
  mids_method <- NULL

  # For every method
  for (m in 1:length(out_cnv[[r]][[1]]$imp_values)) {
    # Is the object null?
    is_null <- is.null(out_cnv[[r]][[1]]$imp_values[[m]])

    if (is_null == FALSE) {
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
            out_cnv$parms$z_m_id,
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
            out_cnv$parms$z_m_id,
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
    } else {
      mids_method[[m]] <- NA
    }
  }

  # Give names to the mids objects
  names(mids_method) <- names(out_cnv[[r]][[1]]$imp_values)

  # Combine repetitions
  mids_r[[r]] <- mids_method
}

# Save slim and prepped mids object
out_cnv <- saveRDS(mids_r, "../output/exp4_conv_20201015_1555_shiny.rds")
