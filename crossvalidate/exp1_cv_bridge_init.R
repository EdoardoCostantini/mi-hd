# Title:    Simulation script for ridge paramter cross-validation
# Project:  Imputing High Dimensional Data
# Author:   Anonymized for peer review
# Created:  2020-08-25
# Modified: 2022-02-23

# Modify parms ------------------------------------------------------------

# Which imputation method to use
parms$meth_sel <- data.frame(DURR_la    = FALSE,
                             IURR_la    = FALSE,
                             blasso     = FALSE,
                             bridge     = TRUE,
                             MI_PCA     = FALSE,
                             MI_CART    = FALSE,
                             MI_RF      = FALSE,
                             MI_qp      = FALSE,
                             MI_am      = FALSE,
                             MI_OP      = FALSE,
                             missFor    = TRUE,
                             mean       = TRUE,
                             CC         = TRUE,
                             GS         = TRUE)

parms$methods <- names(parms$meth_sel)[which(parms$meth_sel==TRUE)]

# What to store (FMI for crossvalidation)
parms$store <-  c(cond         = TRUE,
                  dat_full     = FALSE,
                  dat_miss     = FALSE,
                  sem_EST      = FALSE,
                  sem_CI       = FALSE,
                  lm_EST       = FALSE,
                  lm_CI        = FALSE,
                  fmi          = TRUE,
                  miss_descrps = FALSE,
                  run_time_min = FALSE,
                  imp_values   = FALSE)

# Itereations, repetitions, etc
parms$dt_rep     <- 20
parms$chains     <- 1 
parms$iters      <- 70
parms$burnin_imp <- 50
parms$ndt        <- 10
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
p   <- c(50, 500) # c(50, 500) # number of variables
latent <- c(FALSE, TRUE)[1]
pm <- c(.1, .3)

conds <- expand.grid(ridge, p, latent, pm)#[c(1:2, 10:11, 17:18, 25:26), ]
  colnames(conds) <- c("ridge", "p", "latent", "pm")
