### Title:    Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-05-19

  rm(list = ls())
  library(xtable)
  source("./init_general.R")

## Single Run
  filename <- c("exp1_simOut_20200801_1620", # 750 iterations
                "exp1_simOut_20201130_1006")[2]# 1e3 iterations
  
  # Read R object
  out <- readRDS(paste0("../output/", filename, ".rds"))
  out$parms
  out$conds
  out$session_info
  
## Join multiple results
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
    
  saveRDS(
    output, 
    paste0("../output/", filename, "_res.rds") 
  )
