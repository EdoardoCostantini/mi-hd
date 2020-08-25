### Title:    Initialization scirpt mods for ridge penalty cross validation
### Project:  Imputing High Dimensional Data (exp 2)
### Author:   Edoardo Costantini
### Created:  2020-08-24

# Modify parms ------------------------------------------------------------

# Which imputation method to use
parms$meth_sel <- data.frame(DURR_la = FALSE,
                             DURR_el = FALSE,
                             IURR_la = FALSE,
                             IURR_el = FALSE,
                             bridge  = TRUE,  # the one we want to study
                             blasso  = FALSE,
                             MI_PCA  = FALSE,
                             MI_CART = FALSE,
                             MI_RF   = FALSE,
                             MI_OP   = FALSE,
                             missFor = TRUE,
                             GS      = TRUE,
                             CC      = TRUE)
parms$methods <- names(parms$meth_sel)[which(parms$meth_sel==TRUE)]

# What to store (FMI for crossvalidation)
parms$store <- c(cond         = TRUE,
                 dat_full     = FALSE,
                 dat_miss     = FALSE,
                 sem_EST      = FALSE,
                 sem_CI       = FALSE,
                 lm_EST       = FALSE,
                 lm_CI        = FALSE,
                 fmi          = TRUE, # the thing we are interested in
                 miss_descrps = FALSE,
                 run_time_min = TRUE,
                 imp_values   = FALSE)

# Itereations, repetitions, etc
parms$dt_rep     <- 10
parms$chains     <- 1 
parms$iters      <- 75#75
parms$burnin_imp <- 50#50
parms$ndt        <- 10#10
parms$thin       <- (parms$iters - parms$burnin_imp)/parms$ndt
parms$pos_dt  <- (parms$burnin_imp+1):parms$iters # candidate datasets (after convergence)
parms$keep_dt <- parms$pos_dt[seq(1, 
                                  length(parms$pos_dt), 
                                  parms$thin)] # keep 1 dataset every thin

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
# Modify condtions for crossvalidation

ridge <- 10^seq(from = -1, to = -8, by = -1)

conds <- expand.grid(pm, ridge, r2, int_sub, int_rm, int_da, p)
    
  colnames(conds) <- c("pm", "ridge", "r2", "int_sub", "int_rm", "int_da", "p")
conds <- conds[!(conds$p == p[2] & conds$int_da == TRUE), ]


