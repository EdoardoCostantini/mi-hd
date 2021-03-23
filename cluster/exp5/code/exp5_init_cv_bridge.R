### Title:    Initialization scirpt mods for ridge penalty cross validation
### Project:  Imputing High Dimensional Data (exp 5)
### Author:   Edoardo Costantini
### Created:  2020-07-22


# Baseline ----------------------------------------------------------------
source("./init_general.R")
source("./exp5_init.R")

# Modify parms ------------------------------------------------------------
# Itereations, repetitions, etc
parms$dt_rep     <- 10
parms$chains     <- 1 
parms$iters      <- 2#75
parms$burnin_imp <- 0#50
parms$ndt        <- 2#10
parms$thin       <- (parms$iters - parms$burnin_imp)/parms$ndt
parms$pos_dt  <- (parms$burnin_imp+1):parms$iters # candidate datasets (after convergence)
parms$keep_dt <- parms$pos_dt[seq(1, 
                                  length(parms$pos_dt), 
                                  parms$thin)] # keep 1 dataset every thin

# Always keep missForest to a minimum task
parms$missFor_maxiter <- 2 # maxiter = 20
parms$missFor_ntree <- 10  # ntree = 100

# Of the MI methods, perform only bridge
parms$meth_sel <- list(DURR_all = FALSE,
                       DURR_si  = FALSE,
                       IURR_all = FALSE,
                       IURR_si  = FALSE,
                       blasso   = FALSE,
                       bridge   = TRUE,
                       MI_PCA   = FALSE,
                       MI_CART  = FALSE,
                       MI_RF    = FALSE,
                       MI_OP    = TRUE,
                       missFor  = TRUE,
                       mean     = TRUE,
                       CC       = TRUE,
                       GS       = TRUE)

parms$methods <- names(parms$meth_sel)[which(parms$meth_sel==TRUE)]

# Store only these objects
parms$store <- c(cond     = TRUE,
                 dat_full = FALSE,
                 dat_miss = FALSE,
                 sem_EST  = TRUE,
                 sem_CI   = TRUE,
                 CFA_EST  = TRUE,
                 CFA_CI   = TRUE,
                 semS_EST = TRUE,
                 semS_CI  = TRUE,
                 fmi      = TRUE,
                 miss_des = FALSE,
                 time     = TRUE,
                 time_prep = TRUE,
                 imps     = FALSE)

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
# ridge <- 10^seq(from = -1, to = -8, by = -1) # real values
ridge <- 10^c(-1, -7) # values for trial
    
conds <- expand.grid(ridge, lv, pm, fl, stringsAsFactors = FALSE)
  colnames(conds) <- c("ridge", "lv", "pm", "fl")
