### Title:    Conergence Check Blasso
### Project:  Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-07-03
### Notes:    The goal is to find out how many iterations should be used
###           in the full study for each method.

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
  exp_dat     <- 10
  iters_range <- 1:250 # which set of iterations
  y_range     <- c(2, 3.5)

  # Print ALL
  lapply(seq_along(out_cnv$parms$method[1:8]), function(x){
    mean_traceplot(out_cnv, 
                   method = out_cnv$parms$method[x], 
                   dat = exp_dat, 
                   y_center = TRUE,
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
               y_center = TRUE,
               dat = exp_dat, 
               y_range = y_range, iters = iters_range)

# blasso (200-300 is a good range to chose from)
print(out_cnv$parms$method[3])
iters_range_bl <- 1:250
mean_traceplot(out_cnv,
               method   = out_cnv$parms$method[3],
               dat      = exp_dat, 
               y_range  = c(2, 3.5), iters = iters_range_bl)

# bridge (50-100 good for all)
print(out_cnv$parms$method[4])
mean_traceplot(out_cnv, 
               method   = out_cnv$parms$method[4], 
               dat      = exp_dat,
               y_center = FALSE,
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
               y_center = TRUE,
               dat = exp_dat, 
               y_range = y_range, iters = iters_range)
# MI_RF
print(out_cnv$parms$method[7])
mean_traceplot(out_cnv,
               method = out_cnv$parms$method[7],
               y_center = TRUE,
               dat = exp_dat, 
               y_range = y_range, iters = iters_range)
# MI_T
print(out_cnv$parms$method[8])
mean_traceplot(out_cnv, 
               method = out_cnv$parms$method[8], 
               dat = exp_dat, 
               y_range = y_range, iters = iters_range)

