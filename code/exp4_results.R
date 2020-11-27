### Title:    Results for experiment 4
### Porject:  Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-10-05

  rm(list = ls())
  # library(xtable)
  source("./init_general.R")
  
  # Result selection
  filename <- "exp4_simOut_20201016_2341" # old
  filename <- "exp4_simOut_20201019_1344" # current one
  filename_ad500 <- "exp4_simOut_20201027_1610"
  
  # Read R object
  out <- readRDS(paste0("../output/", filename, ".rds"))
  
  # Describe run
  length(out) - 2
  out$parms$seed
  out$parms$missType
  out$parms$b_int
  out$parms$iters
  
  # or in combine a larger run in two parts
  out_pt1 <- readRDS(paste0("../output/", filename, ".rds"))
  out_pt2 <- readRDS(paste0("../output/", filename_ad500, ".rds"))
  
  # check they were different
  out_pt1[[1]]$`cond_1000_1e-04`$m1_EST$DURR_la
  out_pt2[[1]]$`cond_1000_1e-04`$m1_EST$DURR_la
  
  out_pt1[[400]]$`cond_1000_1e-04`$m1_EST$DURR_la
  out_pt2[[400]]$`cond_1000_1e-04`$m1_EST$DURR_la
  
  # Put together
  out <- c(out_pt1[-c(501:502)], out_pt2)
  
  # fix parms that need to be fixed
  out$parms$dt_rep <- length(out) - 2

# Check presence
  out[[1]]
  names(out[[1]])
  out$conds
  
# Time Analyses -----------------------------------------------------------

  out_time <- sapply(1:length(names(out[[1]])), res_sem_time, out = out)
  colnames(out_time) <- names(out[[1]])
  t(out_time)
  # The warnings are failed convergence methods
  
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

  lapply(which(catch != 8), function(x) out[[x]][[select_cond]]$run_time_min)
  
# Estaimtes Analysis ------------------------------------------------------
  
## Linear Model: Intercept and regression coefficients ##
  m1_res <- lapply(1:length(out[[1]]),
                     function(x) res_sum(out, 
                                         model = "m1", 
                                         condition = x,
                                         bias_sd = TRUE))
  m1_res[[2]]$bias_sd
  par_interest <- c("rel", "trust.s")
  # Show results for a given condition
  m1_bias_per <- lapply(1:length(out[[1]]),
                        function(x) m1_res[[x]]$bias_per)
  m1_bias_esd <- lapply(1:length(out[[1]]),
                        function(x) round(m1_res[[x]]$bias_sd, 1))
  # CI
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
                                       condition = x,
                                       bias_sd = TRUE))
  # Show results for a given condition
  m2_bias_per <- lapply(1:length(out[[1]]),
                        function(x) m2_res[[x]]$bias_per)
  m2_bias_sd <- lapply(1:length(out[[1]]),
                       function(x) round(m2_res[[x]]$bias_sd, 1))
  m2_ci <- lapply(1:length(out[[1]]),
                  function(x) m2_res[[x]]$ci_cov)
  
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
    paste0("../output/", filename_ad500, "_res.rds")
  )
