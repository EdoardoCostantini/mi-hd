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

# DURR --------------------------------------------------------------------
out_DURR <- readRDS("../output/exp5_conv_20210412_1007.rds") # durr
str(out_DURR$parms)
# Methods Names
names(out_DURR[[1]]$`cond_100_0.3_low_1e-07`$imp_values)
methods <- names(out_DURR[[1]]$`cond_100_0.3_low_1e-07`$imp_values)[1:2]

# Iterations and repetations to show
c(Reps   = out_DURR$parms$dt_rep,
  Chians = out_DURR$parms$chains,
  Iters = out_DURR$parms$iters)
exp_dat <- 1 # which data replication
iters_range <- 1:250 # which set of iterations
y_range <- c(4.5, 5.5)

# Actual Plots
mean_traceplot(out_DURR, 
               method = methods[1], 
               dat = 5, 
               y_range = y_range, 
               v_range = 1:10,
               iters = iters_range)

mean_traceplot(out_DURR,
               method = methods[2], 
               dat = 5, 
               y_range = y_range, 
               iters = iters_range)

# IURR_la (50-100 good for all)
mean_traceplot(out_cnv, 
               method = methods[3], 
               dat = exp_dat, 
               y_range = y_range, iters = iters_range)
mean_traceplot(out_cnv, 
               method = methods[4], 
               dat = exp_dat, 
               y_range = y_range, iters = iters_range)

# blasso (try 1000 iterations)
mean_traceplot(out_cnv, 
               method = methods[5], 
               dat = exp_dat, 
               y_range = y_range, iters = iters_range)

# bridge (50-100 good for all)
mean_traceplot(out_cnv, 
               method = methods[6], 
               dat = exp_dat, 
               y_range = y_range, iters = iters_range)


# Tree based MI -----------------------------------------------------------
out_CRF <- readRDS("../output/exp5_conv_20210416_1459.rds") # cart rf
methods <- out_CRF$parms$methods[1:2]
iters_range <- 1:250 # which set of iterations
y_range <- c(4.5, 5.5)

# MI_CART
mean_traceplot(out_CRF, 
               method = methods[1], 
               dat = 3, 
               v_range = 1:10,
               y_range = y_range, iters = iters_range)

# MI_RF
mean_traceplot(out_CRF, 
               method = methods[2], 
               dat = 3,
               v_range = 1:10,
               y_range = y_range, iters = iters_range)

# MI-PCA and MI-OP --------------------------------------------------------
out_PCAOP <- readRDS("../output/exp5_conv_20210416_1457.rds") # mipca
str(out_PCAOP$parms)
methods <- out_PCAOP$parms$methods[1:2]
exp_dat <- 1 # which data replication
iters_range <- 1:250 # which set of iterations
y_range <- c(4.5, 5.5)

# MI_T
mean_traceplot(out_PCAOP, 
               method = methods[1],
               dat = 1,
               v_range = 1:10,
               y_range = y_range, iters = iters_range)

# MI_PCA (50-100 good for all)
mean_traceplot(out_PCAOP, 
               method = methods[2], 
               dat = 5, 
               v_range = 1:3,
               y_range = y_range, 
               iters = iters_range)

