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
exp_dat <- 3 # which data replication (10 possibilities)

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
    dat = exp_dat,
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
