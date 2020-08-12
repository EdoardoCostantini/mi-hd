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
# many chians, for the most difficult condition of experiment 1.

out_cnv <- readRDS("../output/ccheck_exp2_20200812_1115.rds")

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

# What to show
exp_dat <- 3 # which data replication (10 possibilities)
iters_range <- 1:250 # which set of iterations
y_range <- c(-1, .5)

# DURR_la (50-100 good for all)
exp_dat <- 6
mean_traceplot(out_cnv, 
               method = out_cnv$parms$method[1], 
               dat = exp_dat, 
               y_range = y_range, iters = iters_range)
# IURR_la (50-100 good for all)
exp_dat <- 3
mean_traceplot(out_cnv, 
               method = out_cnv$parms$method[2], 
               dat = exp_dat, 
               y_range = y_range, iters = 1:1e3)
# blasso (try 1000 iterations)
exp_dat <- 3
mean_traceplot(out_cnv, 
               method = out_cnv$parms$method[3], 
               dat = exp_dat, 
               y_range = y_range, iters = iters_range)
# Good after 500 for most. Some get better closer to 1e3
# exp_data <- 5 shows very bad z7 and z10, but seem to be 
# exceptions.
# Current proposal: is to perform 1e3 iterations with 950 
# burn in, save 10 data in last 50- exp_dat <- 3 shows why
# we need to get to 1e3 iteartions with variable z5

# bridge (50-100 good for all)
exp_dat <- 3
mean_traceplot(out_cnv, 
               method = out_cnv$parms$method[4], 
               dat = exp_dat,
               y_range = y_range, iters = iters_range)
# MI_PCA (50-100 good for all)
exp_dat <- 1
mean_traceplot(out_cnv,
               method = out_cnv$parms$method[5], 
               dat = exp_dat, 
               y_range = y_range, iters = iters_range)
# MI_CART
exp_dat <- 3
mean_traceplot(out_cnv,
               method = out_cnv$parms$method[6], 
               dat = exp_dat, 
               y_range = y_range, iters = iters_range)
# MI_RF
exp_dat <- 3
mean_traceplot(out_cnv,
               method = out_cnv$parms$method[7],
               dat = exp_dat, 
               y_range = y_range, iters = iters_range)
# MI_T
exp_dat <- 1
mean_traceplot(out_cnv, 
               method = out_cnv$parms$method[8], 
               dat = exp_dat, 
               y_range = y_range, iters = iters_range)

