### Title:    Conergence Check Blasso
### Project:  Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-07-03
### Notes:    The goal is to find out how many iterations should be used
###           in the full study for each method.


# Read Simulation Results -------------------------------------------------
# Special Run with multiple chains. This output does not contain results
# of the simualtion study. It contains an imputation run with many iterations,
# many chians, for the most difficult condition of experiment 1.
out_cnv <- readRDS("../output/cnv_check_20200706_1057.rds")

# What to show
exp_dat <- 5 # which data replication (10 possibilities)
iters_range <- 1:200 # which set of iterations

# DURR_la
mean_traceplot(out_cnv, method = out_cnv$parms$method[1], dat = exp_dat, y_range = c(0,10), iters = iters_range)
# IURR_la
mean_traceplot(out_cnv, method = out_cnv$parms$method[2], dat = exp_dat, y_range = c(0,10), iters = iters_range)
# bridge
mean_traceplot(out_cnv, method = out_cnv$parms$method[3], dat = exp_dat, y_range = c(3,7), iters = iters_range)
# blasso
mean_traceplot(out_cnv, method = out_cnv$parms$method[4], dat = exp_dat, y_range = c(-10,10), iters = iters_range)
mean_traceplot(out_cnv, method = out_cnv$parms$method[4], dat = exp_dat, y_range = c(-2,7), iters = iters_range)
mean_traceplot(out_cnv, method = out_cnv$parms$method[4], dat = exp_dat, y_range = c(-1,1), iters = iters_range)
# MI_PCA
mean_traceplot(out_cnv, method = out_cnv$parms$method[5], dat = exp_dat, y_range = c(0,10), iters = iters_range)
# MI_CART
mean_traceplot(out_cnv, method = out_cnv$parms$method[6], dat = exp_dat, y_range = c(0,10), iters = iters_range)
# MI_RF
mean_traceplot(out_cnv, method = out_cnv$parms$method[7], dat = exp_dat, y_range = c(0,10), iters = iters_range)
# MI_T
mean_traceplot(out_cnv, method = out_cnv$parms$method[8], dat = exp_dat, y_range = c(0,10), iters = iters_range)
