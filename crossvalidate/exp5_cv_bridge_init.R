### Title:    Initialization scirpt mods for ridge penalty cross validation
### Project:  Imputing High Dimensional Data (exp 5)
### Author:   Edoardo Costantini
### Created:  2020-07-22

# Modify parms ------------------------------------------------------------
# Itereations, repetitions, etc
parms$dt_rep     <- 10
parms$chains     <- 1 
parms$iters      <- 75
parms$burnin_imp <- 50
parms$ndt        <- 10
parms$thin       <- (parms$iters - parms$burnin_imp)/parms$ndt
parms$pos_dt  <- (parms$burnin_imp+1):parms$iters # candidate datasets (after convergence)
parms$keep_dt <- parms$pos_dt[seq(1, 
                                  length(parms$pos_dt), 
                                  parms$thin)] # keep 1 dataset every thin

parms$meth_sel <- list(DURR_all = FALSE,   # version w/o SI
                       DURR_si  = FALSE,  # version w/o SI
                       IURR_all = FALSE,   # version w/o SI
                       IURR_si  = FALSE,  # version w/ SI
                       blasso   = TRUE,
                       bridge   = FALSE,
                       MI_PCA   = FALSE,
                       MI_CART  = FALSE,
                       MI_RF    = FALSE,
                       MI_OP    = TRUE,
                       missFor  = TRUE,
                       mean     = TRUE,
                       CC       = TRUE,
                       GS       = TRUE)

parms$methods <- names(parms$meth_sel)[which(parms$meth_sel==TRUE)]

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
lv    <- c(10, 100)       # number of latent variables
pm    <- c(.1, .3)        # proportion of missings level
fl    <- c("high", "low") # factor loadings level
    
conds <- expand.grid(ridge, lv, pm, fl, stringsAsFactors = FALSE)
  colnames(conds) <- c("ridge", "lv", "pm", "fl")
