# Title:    Results for experiment 4
# Project:  Imputing High Dimensional Data
# Author:   Anonymized for peer review
# Created:  2020-10-05
# Modified: 2022-02-16

  rm(list = ls())
  source("./init_general.R")
  
# Single Run -------------------------------------------------------------------

  # filename <- "exp4_simOut_20201204_2121" # updated model 1, 500 data
  # filename <- "exp4_simOut_20201207_1134" # same seed as 20201204_2121, but next 500 samples
  #
  # # Read R object
  # out <- readRDS(paste0("../output/", filename, ".rds"))

# Join additional repetitions multiple results ---------------------------------

  filename1 <- "exp4_simOut_20201204_2121"
  filename2 <- "exp4_simOut_20201207_1134"
  filename <- "exp4_simOut_2020120704" # to save the output
  out_pt1 <- readRDS(paste0("../output/", filename1, ".rds"))
  out_pt2 <- readRDS(paste0("../output/", filename2, ".rds"))
  
  # check they were different
  out_pt1[[1]]$`cond_1000_1e-04`$m1_EST$DURR_la
  out_pt2[[1]]$`cond_1000_1e-04`$m1_EST$DURR_la
  
  out_pt1[[400]]$`cond_1000_1e-04`$m1_EST$DURR_la
  out_pt2[[400]]$`cond_1000_1e-04`$m1_EST$DURR_la
  
  # Put together
  out <- c(out_pt1[-c(501:length(out_pt1))], 
           out_pt2[-c(501:length(out_pt2))])

  # fix parms that need to be fixed
  dt_reps_true <- length(out)
  out$parms <- out_pt1$parms
  out$conds <- out_pt1$conds
  out$parms$dt_rep <- dt_reps_true
  

  # append info from single runs
  out$info <- list(out_pt1 = out_pt1[c(501:length(out_pt1))],
                   out_pt2 = out_pt2[c(501:length(out_pt2))])

# Join additional methods not considered before --------------------------------

# Define file names for the new file to join to the old
nw_filename <- "exp4_simOut_20220131_1603" # new rds filename
filename <- nw_filename

# Read it in of them in R
nw_out <- readRDS(paste0("../output/", nw_filename, ".rds"))

# Extract the meta data from both
meta <- list(og_out = tail(out, 3),
             nw_out = tail(nw_out, 3))

# Get rid of the meta data
og_out <- out[1:out$parms$dt_rep] # temporary mod
nw_out <- nw_out[1:nw_out$parms$dt_rep]

# Append new methods as columns to each repetition and condition

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
rp_filename <- "exp4_simOut_20220226_0950"
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
      # h <- 2
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

# Check it worked
id <- sample(1:1e3, 1)
cbind(og = out[[id]]$`cond_1000_1e-04`$m2_EST[, "bridge"],
      bs = bs_out[[id]]$`cond_1000_1e-04`$m2_EST[, "bridge"],
      rp = rp_out[[id]]$cond_1000_0.01$m2_EST[, "bridge"])

# Fix meta data
out <- bs_out
out <- append(out, meta$rp_out)
out$parms$methods <- unique(c(meta$bs_out$parms$methods, meta$rp_out$parms$methods))


# Time Analyses -----------------------------------------------------------

  out_time <- sapply(1:nrow(out$conds),
                     res_sem_time,
                     out = out,
                     n_reps = out$parms$dt_rep,
                     methods = out$parms$methods[c(1:8, 13:14)]
  )
  colnames(out_time) <- names(out[[1]])
  t(out_time)
  
  # Catch weird runs
  condition <- 2
  select_cond <- names(out[[1]])[condition]
  res_time <- NULL
  catch <- NULL
  
  for (i in 1:out$parms$dt_rep) {
    catch[i] <- length(out[[i]][[select_cond]]$run_time_min)
    res_time <- rbind(res_time, out[[i]][[select_cond]]$run_time_min)
  }
  which(catch != 8)

  # Detect Defective Method
  # Regular Runs
  lapply(which(catch == 8), function(x) out[[x]][[select_cond]]$run_time_min)[1]
  # Weird Runs
  lapply(which(catch != 8), function(x) out[[x]][[select_cond]]$run_time_min)
  
# Estimates Analysis ------------------------------------------------------
  
## Linear Model: Intercept and regression coefficients ##
  m1_res <- lapply(1:length(out[[1]]),
                     function(x) res_sum(out, 
                                         model = "m1", 
                                         condition = x,
                                         bias_sd = TRUE))

  # Bias for a given model
  lapply(1:length(out[[1]]), function(x) m1_res[[x]]$bias_per)
  lapply(1:length(out[[1]]), function(x) round(m1_res[[x]]$bias_sd, 1))
  
  # CI
  lapply(1:length(out[[1]]), function(x) m1_res[[x]]$ci_cov)
  
  # Available results for each method (converged MI algorithm)
  m1_res[[1]]$validReps
  m1_res[[2]]$validReps
  
## Linear Model 2 ##
  m2_res <- lapply(1:length(out[[1]]),
                   function(x) res_sum(out, 
                                       model = "m2", 
                                       condition = x,
                                       bias_sd = TRUE))
  # Bias for model 2
  lapply(1:length(out[[1]]), function(x) m2_res[[x]]$bias_per)
  lapply(1:length(out[[1]]), function(x) round(m2_res[[x]]$bias_sd, 1))
  
  # CI
  lapply(1:length(out[[1]]), function(x) m2_res[[x]]$ci_cov)
  
  # Available results
  m2_res[[1]]$validReps
  m2_res[[2]]$validReps
  
## Overall Euclidean Distances
  ed_all_res <- lapply(1:length(out[[1]]),
                       function(x) res_ed_overall(out,
                                                  condition = x))
  
# Save Results ------------------------------------------------------------
  output <- lapply(list(m1 = m1_res,
                        m2 = m2_res,
                        ed_all = ed_all_res), 
                   function(x){
                     names(x) <- paste0("cond", seq_along(out[[1]]))
                     return(x)
                   }
  )
  output$parms <- out$parms
  output$conds <- out$conds
  
  saveRDS(
    output, 
    paste0("../output/", "exp4_simOut_20220226_0950", "_res.rds")
  )
