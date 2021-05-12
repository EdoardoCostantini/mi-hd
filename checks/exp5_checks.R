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
  
# Factor Loading Size on EVS ----------------------------------------------

  library(semPlot)
  library(lavaan)
  
  rm(list=ls())
  source("./init_general.R")
  source("./exp5_init.R")
  
  # Read EVS data
  dat <- readRDS("../data/exp4_EVS2017_full.rds")$full
  
  # Define CFA model (choose 1 V)
  # V1
  item_ids <- list(trustState = c(120, 121, 127, 131))
  # V2
  item_ids <- list(familyVals = 82:84,
                   genderRole = 72:79)
  # V3
  item_ids <- list(poliAction = 98:101,
                   natatt = 185:187)
  
  # For the chosen group of variables, do the following
  lv_names <- names(item_ids)
  lv_number = length(item_ids)
  
  lv_models <- sapply(1:lv_number, function(i){
    paste0(lv_names[i], 
           " =~ ",
           paste0("v", item_ids[[i]], collapse = " + ")
    )
    
  })
  CFA_model <- paste(lv_models, collapse = "\n")
  
  # Fit CFA
  fit <- cfa(CFA_model, data = dat, std.lv = TRUE)
  
  # Evaluate Fit
  fit.out <- summary(fit, fit.measures = TRUE, standardized = TRUE)
  round(fit.out$FIT[c("cfi", "tli", "rmsea", "rmsea.pvalue")], 3)
  
  # Compare
  CFA_par <- parameterEstimates(fit, 
                                se = FALSE, zstat = FALSE, 
                                pvalue = FALSE, ci = FALSE,
                                standardized = TRUE)
  
  # Measured
  # Factor loadings
  CFA_par[1:length(unlist(item_ids)), c(1:3, 6)]
  
  # Plot
  semPaths(fit, "std",
           nCharNodes = max(nchar(names(item_ids)))
           )
           