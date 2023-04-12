# Project:   imputeHD-comp
# Objective: Process result files to prepare analysis
# Author:    Edoardo Costantini
# Created:   2023-04-12
# Modified:  2023-04-12
# Notes: 

rm(list = ls())
source("./init_general.R")

# Load data --------------------------------------------------------------------

filename <- c(
    "exp1_2_simOut_20230408_1748" # 30 reps
)[1]

# Read R object
out <- readRDS(paste0("../output/", filename, ".rds"))
out$parms
out$conds
out$session_info

# Time Analyses -----------------------------------------------------------

  out_time <- sapply(1:nrow(out$conds),
                     res_sem_time,
                     out = out,
                     n_reps = out$parms$dt_rep,
                     methods = out$parms$methods[c(1:11)]
  )

  colnames(out_time) <- names(out[[1]])
  t(out_time)

  # Save
  saveRDS(
      t(out_time),
      paste0("../output/", filename, "_time.rds")
  )
  
# Univariate Analyses -----------------------------------------------------
  
## MLE estimates (saturated sem model) ##

  # Extract results per conditions
  sem_res <- lapply(1:length(out[[1]]),
                    function(x) res_sum(out, 
                                        model = "sem",
                                        condition = x))

  # Show results for a given condition
  lapply(1:length(out[[1]]),
         function(x) sem_res[[x]]$bias_per)

## Linear Model: Intercept and regression coefficients ##

  lm_res <- lapply(1:length(out[[1]]),
                   function(x) res_sum(out, 
                                       model = "lm", 
                                       condition = x))
  
  # Show results for a given condition
  lapply(1:length(out[[1]]),
         function(x) lm_res[[x]]$bias_per)

# Save Results ------------------------------------------------------------

    output <- lapply(
        list(
            sem = sem_res,
            lm = lm_res
        ),
        function(x) {
            names(x) <- paste0("cond", seq_along(out[[1]]))
            return(x)
        }
    )
    output$parms <- out$parms
    output$conds <- out$conds

    # Transform for plot lm model
    gg_out_lm <- plotwise(
        res = output,
        model = "lm",
        parPlot = list(Betas = 1:5),
        item_group = c(1:3), # items in a group receiving miss values
        meth_compare = c(
            "DURR_la", "IURR_la", "blasso", "bridge",
            "MI_PCA", "MI_CART", "MI_RF", "stepFor",
            "MI_OP",
            "CC", "GS", "MI_qp", "MI_am"
        ),
        exp_factors = c("p", "collinearity")
    )

    # Save
    saveRDS(
        gg_out_lm,
        paste0("../output/", filename, "_res_lm.rds")
    )

    # Transform for plot sem model
    gg_out_sem <- plotwise(
        res = output,
        model = "sem",
        parPlot = list(
            Means = 1:6,
            Variances = 7:12,
            Covariances = 13:27
        ),
        item_group = c(1:3), # items in a group receiving miss values
        meth_compare = c(
            "DURR_la", "IURR_la", "blasso", "bridge",
            "MI_PCA", "MI_CART", "MI_RF", "stepFor",
            "MI_OP",
            "CC", "GS", "MI_qp", "MI_am"
        ),
        exp_factors = c("p", "collinearity")
    )

    # Save
    saveRDS(
        gg_out_sem,
        paste0("../output/", filename, "_res.rds")
    )