### Title:    Results for experiment 4
### Project:  Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-10-05

  rm(list = ls())
  source("./init_general.R")
  
# Analysis single result
  filename <- "exp4_simOut_20201204_2121" # updated model 1, 500 data
  filename <- "exp4_simOut_20201207_1134" # same seed as 20201204_2121, but next 500 samples
  
  # Read R object
  out <- readRDS(paste0("../output/", filename, ".rds"))

## Join multiple results
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
  
# Time Analyses -----------------------------------------------------------

  out_time <- sapply(1:length(names(out[[1]])), res_sem_time, out = out)
  colnames(out_time) <- names(out[[1]])
  t(out_time)
  # The warnings are failed convergence methods
  
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
    paste0("../output/", filename, "_res.rds")
  )
