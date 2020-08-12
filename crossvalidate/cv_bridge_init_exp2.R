### Title:    Initialization scirpt mods for ridge penalty cross validation
### Project:  Imputing High Dimensional Data (exp 2)
### Author:   Edoardo Costantini
### Created:  2020-07-22

# Modify parms ------------------------------------------------------------
# Itereations, repetitions, etc
parms$dt_rep     <- 10  # reduced number of data repretions just to make sure 
# that ridge penalty decision is not done on 1 dataset

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
                             CC      = TRUE
)
parms$methods <- names(parms$meth_sel)[which(parms$meth_sel==TRUE)]

# Report names
parms$report_file_name <- paste0("cv_bridge_",
                                 "exp", parms$exp, "_",
                                 parms$start_time, 
                                 ".txt")
parms$results_file_name <- paste0("cv_bridge_",
                                  "exp", parms$exp, "_",
                                  parms$start_time,
                                  ".rds")
parms$description <- c("In each repetition, 1 dataset is created for each condition.
        Imputation methods are used on that condition-specific dataset.
        Results are therefore given per dataset in condition")


# Conditions --------------------------------------------------------------
# Modify condtions for crossvalidation
ridge <- c(1e-2, 1e-3, 1e-4, 1e-5, 1e-6)
lv    <- c(10, 100)       # number of latent variables
pm    <- c(.1, .3)        # proportion of missings level
fl    <- c("high", "low") # factor loadings level
    
conds <- expand.grid(ridge, lv, pm, fl, stringsAsFactors = FALSE)
  colnames(conds) <- c("ridge", "lv", "pm", "fl")
