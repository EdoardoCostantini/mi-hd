# Title:    Putting result object together
# Project:  Imputing High Dimensional Data
# Author:   Edoardo Costantini
# Created:  2020-05-19
# Modified: 2022-02-04

  rm(list = ls())
  source("./init_general.R")

# Single Run -------------------------------------------------------------------
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
nw_filename <- "exp1_simOut_20220201_1749" # new rds filename

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

# Replace results of a method that was re-run ----------------------------------

# Load results
rp_filename <- "exp1_simOut_20220214_1418"
rp_out <- readRDS(paste0("../output/", rp_filename, ".rds"))

# Extract the meta data from both
meta <- list(bs_out = tail(out, 3), # base file
             rp_out = tail(rp_out, 3)) # replacing file

# Get rid of the meta data
bs_out <- out[1:out$parms$dt_rep] # temporary mod
rp_out <- rp_out[1:rp_out$parms$dt_rep]

# Original order of methods
mt_order <- colnames(bs_out[[1]]$cond_50_FALSE_0.1$sem_EST)
mt_order_mitime <- names(bs_out[[1]]$cond_50_FALSE_0.1$run_time_min)

# Append new methods as columns to each repetition and condition --------
rp_method <- c("bridge")

for(i in 1:length(bs_out)){ # for every repetition
  # i <- 1
  for(j in 1:length(bs_out[[i]])){ # for every condition
    # j <- 1
    for(h in 2:length(bs_out[[i]][[j]])){ # for every object stored (excpet condition label)
      # h <- 6
      print(paste0("i = ", i, "; j = ", j, "; h = ", h))
      multi_dim <- length(dim(bs_out[[i]][[j]][[h]])) == 2
      if(multi_dim){
        # Check method was succesful
        bs_success <- rp_method %in% colnames(bs_out[[i]][[j]][[h]])
        rp_success <- rp_method %in% colnames(rp_out[[i]][[j]][[h]])
        if(bs_success & rp_success){
          # Check method was successful in run
          bs_out[[i]][[j]][[h]][, rp_method] <-
            rp_out[[i]][[j]][[h]][, rp_method]
        }
        if(!bs_success & rp_success){
          print(paste0("i = ", i, "; j = ", j, "; h = ", h))
          # Append column
          bs_out[[i]][[j]][[h]] <- cbind(bs_out[[i]][[j]][[h]],
                                         rp_out[[i]][[j]][[h]][, rp_method, drop = FALSE])
          # Put in original order
          bs_out[[i]][[j]][[h]] <- bs_out[[i]][[j]][[h]][, mt_order]
        }
        if(bs_success & !rp_success){
          colindex <- !colnames(bs_out[[i]][[j]][[h]]) %in% rp_method
          bs_out[[i]][[j]][[h]] <- bs_out[[i]][[j]][[h]][, colindex]
        }
        if(!bs_success & !rp_success){
          next
        }
      } else {
        # Check method was succesful
        bs_success <- rp_method %in% names(bs_out[[i]][[j]][[h]])
        rp_success <- rp_method %in% names(rp_out[[i]][[j]][[h]])
        if(bs_success & rp_success){
          # Replcae column
          bs_out[[i]][[j]][[h]][rp_method] <-
            rp_out[[i]][[j]][[h]][rp_method]
        }
        if(!bs_success & rp_success){
          print(paste0("i = ", i, "; j = ", j, "; h = ", h))
          # Append value
          bs_out[[i]][[j]][[h]] <- c(bs_out[[i]][[j]][[h]],
                                     rp_out[[i]][[j]][[h]][rp_method])
          # Put in original order
          bs_out[[i]][[j]][[h]] <- bs_out[[i]][[j]][[h]][mt_order_mitime]
        }
        if(bs_success & !rp_success){
          # Get rid of column
          index <- !names(bs_out[[i]][[j]][[h]]) %in% rp_method
          bs_out[[i]][[j]][[h]] <- bs_out[[i]][[j]][[h]][index]
        }
        if(!bs_success & !rp_success){
          next
        }
      }
    }
  }
}

out <- bs_out
out <- append(out, meta$rp_out)
out$parms$methods <- unique(c(meta$bs_out$parms$methods, meta$rp_out$parms$methods))

id <- sample(1:1e3, 1)
cbind(bs_out[[id]]$cond_50_FALSE_0.3$sem_EST[, "bridge"],
      rp_out[[id]]$cond_50_FALSE_0.3$sem_EST[, "bridge"])

# Time Analyses -----------------------------------------------------------

  out_time <- sapply(1:nrow(out$conds),
                     res_sem_time,
                     out = out,
                     n_reps = out$parms$dt_rep,
                     methods = out$parms$methods[c(1:8, 13:14)]
  )

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
  filename <- paste0("exp1_simOut_20220214_1418")
  saveRDS(
    gg_out_sem,
    paste0("../output/", filename, "_res.rds")
  )