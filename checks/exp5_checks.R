### Title:    Checking parts of experiment 5 (addendum)
### Project:  Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2021-01-29

# Data Generation Latent Structure ----------------------------------------

  rm(list=ls())
  source("./init_general.R")
  source("./exp2_init.R")

#  > Bridge is faster without robustness checks ---------------------------

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

