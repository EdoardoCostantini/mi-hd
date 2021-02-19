### Title:    Checking parts of experiment 5 (addendum)
### Project:  Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2021-01-29

# Bridge is faster without robustness checks ------------------------------

  rm(list=ls())
  source("./init_general.R")
  source("./exp2_init.R")

  cond <-  conds[2, ]

## Data ------------------------------------------------------------------ ##
  
  set.seed(20210127)
  
  simData_list <- simData_lv(parms, cond)
  Xy     <- simData_list$dat
  Xy_mis <- imposeMiss_lv_MD(simData_list, parms, cond, plot = TRUE)
  O <- data.frame(!is.na(Xy_mis)) # matrix index of observed values

# Robustness checks ON
  bridge <- impute_BRIDGE(Z = Xy_mis,
                              O = O,
                              ridge_p = cond$ridge,
                              parms = parms,
                              perform = parms$meth_sel$bridge,
                              robust = TRUE)
  # When dimensionality is high, Bridge often end in singularity issue
  # After 1 or two iterations
  bridge$time

# Robustness checks OFF
  bridge.fast <- impute_BRIDGE(Z = Xy_mis,
                               O = O,
                               ridge_p = cond$ridge,
                               parms = parms,
                               perform = parms$meth_sel$bridge,
                               robust = FALSE)
  bridge.fast$time
  # As there is no need for these checks with continuos data
  # then I'll use this version because it's three times faster!

# The method will also run with preprocessing
  prep_Si_out <- prep_SI(dt_in = Xy_mis,
                         model.var = parms$z_m_id,
                         maxit = 10)
  Xy_SI <- prep_Si_out$dt_out
  
  bridge_PRE <- impute_BRIDGE(Z = Xy_SI,
                              O = data.frame(!is.na(Xy_SI)),
                              ridge_p = cond$ridge,
                              parms = parms,
                              perform = parms$meth_sel$bridge,
                              robust = TRUE)
  bridge_PRE$time

# Computation Time --------------------------------------------------------

  out <- readRDS("../output/exp5_simOut_20210201_1850.rds")
  
  iters_test <- out$parms$iters
  iters_goal <- 60
    
  sapply(1:out$parms$dt_rep, function(cond_i){
    # cond_i <- 1
    imp_time <- t(sapply(1:nrow(out$conds), function(i){
      out[[cond_i]][[i]]$run_time_min
    }))
    
    prp_time <- sapply(1:nrow(out$conds), function(i){
      out[[cond_i]][[i]]$run_time_prep # not exactly check what it is
    })
    
    sum( rowSums(imp_time * iters_goal/iters_test) ) / 60
  })
  
  # Estimate for LISA Cluster
  # Lisa Cluster computes cpu hours as: 
  # number of nodes *  number of cores *  job time in hours.
  n_nodes <- 4
  n_cores <- 16 # fixed ?
  job_time <- 25 # hours
  n_nodes * n_cores * job_time
  
  # I need to perform 1000 dt reps
  1e3
  # And I can perform 15 of them (16 cores - 1) per node
  1e3 / 15
  # So I need about 67 nodes
  
  n_nodes <- 70
  n_cores <- 16 # fixed ?
  job_time <- 25 # hours
  n_nodes * n_cores * job_time
  