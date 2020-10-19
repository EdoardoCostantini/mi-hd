### Title:    Initialization scirpt mods for ridge penalty cross validation
### Project:  Imputing High Dimensional Data (exp 2)
### Author:   Edoardo Costantini
### Created:  2020-08-24

# Modify parms ------------------------------------------------------------

# Which imputation method to use
parms$meth_sel <- data.frame(DURR_la    = FALSE,
                             IURR_la    = FALSE,
                             blasso     = FALSE,
                             bridge     = TRUE,
                             MI_PCA     = FALSE,
                             MI_CART    = FALSE,
                             MI_RF      = FALSE,
                             MI_OP      = FALSE,
                             missFor    = TRUE,
                             mean       = TRUE,
                             GS         = TRUE,
                             CC         = TRUE)

parms$methods <- names(parms$meth_sel)[which(parms$meth_sel==TRUE)]

# What to store (FMI for crossvalidation)
parms$store <- c(cond         = TRUE,
                 dat_full     = FALSE,
                 dat_miss     = FALSE,
                 m1_EST       = FALSE,
                 m1_CI        = FALSE,
                 m2_EST       = FALSE,
                 m2_CI        = FALSE,
                 fmi          = TRUE, # the thing we are interested in
                 miss_descrps = FALSE,
                 run_time_min = TRUE,
                 imp_values   = FALSE)

# Itereations, repetitions, etc
parms$dt_rep     <- 3 # 20
parms$chains     <- 1 
parms$iters      <- 5 # 1e2
parms$burnin_imp <- 0 # 60
parms$ndt        <- 5 # 50
parms$thin       <- (parms$iters - parms$burnin_imp)/parms$ndt
parms$pos_dt     <- (parms$burnin_imp+1):parms$iters
parms$keep_dt    <- parms$pos_dt[seq(1, 
                                     length(parms$pos_dt), 
                                     parms$thin)]

# Report names
parms$report_file_name <- paste0("exp",
                                 parms$exp, "_",
                                 "cv_bridge_",
                                 parms$start_time,
                                 ".txt")
parms$results_file_name <- paste0("exp",
                                  parms$exp, "_",
                                  "cv_bridge_",
                                  parms$start_time,
                                  ".rds")
parms$description <- c("In each repetition, 1 dataset is created for each condition.
        Imputation methods are used on that condition-specific dataset.
        Results are therefore given per dataset in condition")


# Conditions --------------------------------------------------------------

  # Fixed random factor
  parms$pm <- c(.1, .2) # make this a random factor
  
  # Experimental factors
  n     <- c(1e3, 3e2) # number of observations
  
  # Dataframe of conditions w/ ridge
  ridge <- 10^seq(from = 4, to = -8, by = -1)
  conds <- expand.grid(n = n, ridge = ridge)
  
