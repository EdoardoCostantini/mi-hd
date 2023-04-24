# Project:   imputeHD-comp
# Objective: Process result files to prepare analysis
# Author:    Edoardo Costantini
# Created:   2023-04-12
# Modified:  2023-04-21
# Notes: 

rm(list = ls())
source("./init_general.R")

# Load R methods except MI-QP results ------------------------------------------

filename <- c(
    "exp1_2_simOut_20230408_1748", # 30 reps with MI-QP time estimate
    "exp1_2_simOut_20230419_1403" # All R methods expect MI-QP with correct collinearity
)[2]

# Read R object
out_Rmeth <- readRDS(paste0("../output/", filename, ".rds"))

# Load IVEware results ---------------------------------------------------------

filename <- c(
    "exp1_2_simOut_20230421_1151"
)[1]

# Read R object
out_IVEware <- readRDS(paste0("../output/", filename, ".rds"))

# Load MI-QP results -----------------------------------------------------------

filename <- c(
    "exp1_2_simOut_20230419_0957"
)[1]

# Read R object
out_MIQP <- readRDS(paste0("../output/", filename, ".rds"))

# Final file name --------------------------------------------------------------

filename <- paste0("exp1_2_simOut_", format(Sys.time(), "%Y%m%d_%H%M"))

# Concatenate results ----------------------------------------------------------

# Add IVEware to Rmeth

for (i in 1:(out_Rmeth$parms$dt_rep)) { # for every repetition
    # i <- 1
    for (j in 1:(nrow(out_Rmeth$conds))) { # for every condition
        # j <- 1
        for (h in 2:length(out_Rmeth[[i]][[j]])) {
            # h <- 2
            multi_dim <- length(dim(out_Rmeth[[i]][[j]][[h]])) == 2
            if (multi_dim) {
                colnames1 <- colnames(out_Rmeth[[i]][[j]][[h]])
                colnames2 <- colnames(out_IVEware[[i]][[j]][[h]])
                colindex <- !colnames2 %in% colnames1
                out_Rmeth[[i]][[j]][[h]] <- cbind(
                    out_Rmeth[[i]][[j]][[h]],
                    out_IVEware[[i]][[j]][[h]][, colindex,
                        drop = FALSE
                    ]
                )
            } else {
                names1 <- names(out_Rmeth[[i]][[j]][[h]])
                names2 <- names(out_IVEware[[i]][[j]][[h]])
                namesindex <- !names2 %in% names1
                out_Rmeth[[i]][[j]][[h]] <- c(
                    out_Rmeth[[i]][[j]][[h]],
                    out_IVEware[[i]][[j]][[h]][namesindex]
                )
            }
        }
    }
}

# Add MI-QP

for (i in 1:(out_Rmeth$parms$dt_rep)) { # for every repetition
    # i <- 1
    for (j in 1:nrow(out_MIQP$conds)) { # for every condition
        # j <- 4
        for (h in 2:length(out_Rmeth[[i]][[j]])) {
            # h <- 2
            multi_dim <- length(dim(out_Rmeth[[i]][[j]][[h]])) == 2
            if (multi_dim) {
                colnames1 <- colnames(out_Rmeth[[i]][[j]][[h]])
                colnames2 <- colnames(out_MIQP[[i]][[j]][[h]])
                colindex <- !colnames2 %in% colnames1
                out_Rmeth[[i]][[j]][[h]] <- cbind(
                    out_Rmeth[[i]][[j]][[h]],
                    out_MIQP[[i]][[j]][[h]][, colindex,
                        drop = FALSE
                    ]
                )
            } else {
                names1 <- names(out_Rmeth[[i]][[j]][[h]])
                names2 <- names(out_MIQP[[i]][[j]][[h]])
                namesindex <- !names2 %in% names1
                out_Rmeth[[i]][[j]][[h]] <- c(
                    out_Rmeth[[i]][[j]][[h]],
                    out_MIQP[[i]][[j]][[h]][namesindex]
                )
            }
        }
    }
}

# Update the parms objects -----------------------------------------------------

# Start from the first file's parms object
parms <- out_Rmeth$parms

# Create new vector of methds
new_meth <- rep(TRUE, length(parms$meth_sel))

# Give it good names
names(new_meth) <- names(parms$meth_sel)

# Update internal objects 
parms$meth_sel <- new_meth
parms$methods <- names(new_meth)

# Create final output object
out <- out_Rmeth
out$parms <- parms

# Time Analyses ----------------------------------------------------------------

  out_time <- sapply(1:nrow(out$conds),
                     res_sem_time,
                     out = out,
                     n_reps = out$parms$dt_rep,
                     methods = out$parms$methods#[c(1:11)]
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
            "DURR_la",
            "IURR_la",
            "blasso",
            "bridge",
            "MI_PCA",
            "MI_CART",
            "MI_RF",
            "stepFor",
            "MI_OP",
            "CC",
            "GS",
            "MI_qp",
            "MI_am"
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
            "DURR_la", 
            "IURR_la", 
            "blasso", 
            "bridge",
            "MI_PCA",
            "MI_CART",
            "MI_RF",
            "stepFor",
            "MI_OP",
            "CC",
            "GS",
            "MI_qp",
            "MI_am"
        ),
        exp_factors = c("p", "collinearity")
    )

    # Save
    saveRDS(
        gg_out_sem,
        paste0("../output/", filename, "_res.rds")
    )