### Title:    Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-05-19

  rm(list = ls())
  library(xtable)
  source("./init_general.R")
  
  # Checks
  # filename <- "sim_res_20200801_1620" # the 750 data reps and 100 iters run
  # 20200715 Correct
  # filename <- "sim_res_20200710_1019" # first submision
  # Current Correct
  
  # "exp1_simOut_20200801_1620"
  filename <- c("sim_res_20200710_1019",
                "exp1_simOut_20200731_1735", # <- CURRENT
                "exp1_simOut_20200801_1620")[3]
  # Result 1 is the one used for first submission to EAM
  # Result 2 and 3 corrected for MNAR issue. They are equivalent, but
  # second file has more repetitions (500 vs 750)
  
  # Read R object
  out <- readRDS(paste0("../output/", filename, ".rds"))

# Time Analyses -----------------------------------------------------------

  out_time <- sapply(1:length(names(out[[1]])), res_sem_time, out = out)
  colnames(out_time) <- names(out[[1]])
  t(out_time)
  
# Univariate Analyses -----------------------------------------------------
  
## MLE estiamtes (saturated sem model) ##
  
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
  output$conds <- conds
    
  saveRDS(
    output, 
    paste0("../output/", filename, "_res.rds") 
  )
