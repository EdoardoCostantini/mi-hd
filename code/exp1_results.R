# Title:    Putting result object together
# Project:  Imputing High Dimensional Data
# Author:   Edoardo Costantini
# Created:  2020-05-19
# Modified: 2022-01-31

  rm(list = ls())
  source("./init_general.R")

## Single Run
  filename <- c("exp1_simOut_20200801_1620", # 750 iterations
                "exp1_simOut_20201130_1006", # 1e3 iterations
                "exp1_simOut_20220128_1247")[3]
  
  # Read R object
  out <- readRDS(paste0("../output/", filename, ".rds"))
  out$parms
  out$conds
  out$session_info

# Join additional repetitions multiple results ----------------------------

  filename1 <- "exp1_simOut_20201215_1018"
  filename2 <- "exp1_simOut_20201216_1645"
  filename <- "exp1_simOut_2020121516" # to save the output
  
  # Read R objects
  out_pt1 <- readRDS(paste0("../output/", filename1, ".rds"))
  out_pt2 <- readRDS(paste0("../output/", filename2, ".rds"))
  
  # Put together
  out <- c(out_pt1[-c((out_pt1$parms$dt_rep+1):length(out_pt1))], 
           out_pt2[-c((out_pt1$parms$dt_rep+1):length(out_pt2))])
  
  # fix parms that need to be fixed
  dt_reps_true <- length(out)
  out$parms <- out_pt1$parms
  out$conds <- out_pt1$conds
  out$parms$dt_rep <- dt_reps_true
  
  # append info from single runs
  out$info <- list(out_pt1 = out_pt1[c(501:length(out_pt1))],
                   out_pt2 = out_pt2[c(501:length(out_pt2))])

# Join additional methods not considered before

# Define file names for the two files to join
og_filename <- "exp1_simOut_20201130_1006" # original rds filename
nw_filename <- "exp1_simOut_20220128_1635" # new rds filename

# Read both of them in R
og_out <- readRDS(paste0("../output/", og_filename, ".rds"))
nw_out <- readRDS(paste0("../output/", nw_filename, ".rds"))

# Extract the meta data from both
meta <- list(og_out = tail(og_out, 3),
             nw_out = tail(nw_out, 3))

# Get rid of the meta data
og_out <- og_out[1:og_out$parms$dt_rep] # temporary mod
nw_out <- nw_out[1:nw_out$parms$dt_rep]

# Append new methods as columns to each repetition and condition --------

for(i in 1:length(og_out)){ # for every repetition
  for(j in 1:length(og_out[[i]])){ # for every condition
    for(h in 2:length(og_out[[i]][[j]])){
      multi_dim <- length(dim(og_out[[i]][[j]][[h]])) == 2
      if(multi_dim){
        colnames1 <- colnames(og_out[[i]][[j]][[h]])
        colnames2 <- colnames(nw_out[[i]][[j]][[h]])
        colindex <- !colnames2 %in% colnames1
        og_out[[i]][[j]][[h]] <- cbind(og_out[[i]][[j]][[h]],
                                       nw_out[[i]][[j]][[h]][, colindex,
                                                               drop = FALSE])
      } else {
        names1 <- names(og_out[[i]][[j]][[h]])
        names2 <- names(nw_out[[i]][[j]][[h]])
        namesindex <- !names2 %in% names1
        og_out[[i]][[j]][[h]] <- c(og_out[[i]][[j]][[h]],
                                   nw_out[[i]][[j]][[h]][namesindex])
      }
    }
  }
}

out <- og_out
out <- append(out, meta$nw_out)
out$parms$methods <- unique(c(meta$og_out$parms$methods, meta$nw_out$parms$methods))

# Time Analyses -----------------------------------------------------------

  out_time <- sapply(1:length(names(out[[1]])), res_sem_time, out = out)
  colnames(out_time) <- names(out[[1]])
  t(out_time)
  
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

  output <- lapply(list(sem   = sem_res,
                        lm    = lm_res), 
                   function(x){
                     names(x) <- paste0("cond", seq_along(out[[1]]))
                     return(x)
                   }
  )
  output$parms <- out$parms
  output$conds <- out$conds

  # Transform for plot
  gg_out_sem <- plotwise(res = output,
                         model = "sem",
                         parPlot = list(Means = 1:6,
                                        Variances = 7:12,
                                        Covariances = 13:27),
                         item_group = c(1:3), # items in a group recieving miss values
                         meth_compare = c("DURR_la","IURR_la", "blasso", "bridge",
                                          "MI_PCA", "MI_CART" ,"MI_RF",
                                          "MI_OP",
                                          "CC", "GS", "MI_qp", "MI_am"))

  # Save
  filename <- paste0("exp1_simOut_20220128_1635_", "joined_", "20201130_1006")
  saveRDS(
    gg_out_sem,
    paste0("../output/", filename, "_res.rds")
  )