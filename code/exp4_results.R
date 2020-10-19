### Title:    Results for experiment 4
### Porject:  Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-10-05

  rm(list = ls())
  library(xtable)
  source("./init_general.R")
  
  # Result selection
  filename <- "exp4_simOut_20201013_1846" # a full run with an issue at data 2
  filename <- "exp4_simOut_20201015_1227" # a full run with few data repetitions
  filename <- "exp4_simOut_20201015_1553"
  filename <- "exp4_simOut_20201016_2341" # current one
  
  # Read R object
  out <- readRDS(paste0("../output/", filename, ".rds"))
  
  # Describe run
  out$parms$missType
  out$parms$b_int
  out$parms$iters
  
# Check presence
  out[[1]]
  names(out[[1]])
  out$conds
  
# Time Analyses -----------------------------------------------------------

  out_time <- sapply(1:length(names(out[[1]])), res_sem_time, out = out)
  colnames(out_time) <- names(out[[1]])
  t(out_time)
  
  # Catch wierd runs
  # Time
  condition <- 2
  select_cond <- names(out[[1]])[condition]
  res_time <- NULL
  catch <- NULL
  
  for (i in 1:out$parms$dt_rep) {
    catch[i] <- length(out[[i]][[select_cond]]$run_time_min)
    res_time <- rbind(res_time, out[[i]][[select_cond]]$run_time_min)
  }
  which(catch != 8)
  length(which(catch != 8))
  
  out[[21]][[select_cond]]$run_time_min
  out[[91]][[select_cond]]$run_time_min
  out[[105]][[select_cond]]$run_time_min
  out[[429]][[select_cond]]$run_time_min
  out[[428]][[select_cond]]$run_time_min
  
# Estaimtes Analysis ------------------------------------------------------
  
## Linear Model: Intercept and regression coefficients ##
  m1_res <- lapply(1:length(out[[1]]),
                     function(x) res_sum(out, 
                                         model = "m1", 
                                         condition = x))
  
  # Show results for a given condition
  m1_bias <- lapply(1:length(out[[1]]),
                    function(x) m1_res[[x]]$bias_per)
  m1_ci <- lapply(1:length(out[[1]]),
                  function(x) m1_res[[x]]$ci_cov)
  
  # Available results for each method (converged MI algorithm)
  m1_res[[1]]$validReps
  m1_res[[2]]$validReps
    # 10 MI_OP and 1 bridge failed
  
## Linear Model 2 ##
  m2_res <- lapply(1:length(out[[1]]),
                   function(x) res_sum(out, 
                                       model = "m2", 
                                       condition = x))
  # Show results for a given condition
  m2_bias <- lapply(1:length(out[[1]]),
                    function(x) m2_res[[x]]$bias_per)
  m2_ci <- lapply(1:length(out[[1]]),
                  function(x) m2_res[[x]]$ci_cov)
  
  # Available results
  m2_res[[1]]$validReps
  m2_res[[2]]$validReps

# Save Results ------------------------------------------------------------
  output <- lapply(list(m1 = m1_res,
                        m2 = m2_res), 
                   function(x){
                     names(x) <- paste0("cond", seq_along(out[[1]]))
                     return(x)
                   }
  )
  
  saveRDS(
    output, 
    paste0("../output/", filename, "_res.rds")
  )
