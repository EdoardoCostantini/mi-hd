### Title:    Results for experiment 3
### Porject:  Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-08-25

  rm(list = ls())
  library(xtable)
  source("./init_general.R")
  
  # Result selection
  filename <- "exp3_simOut_20200825_0951"
  
  # Read R object
  out <- readRDS(paste0("../output/", filename, ".rds"))
  
# Check presence
  out[[1]]
  names(out[[1]])
  out$conds
  
# Time Analyses -----------------------------------------------------------

  out_time <- sapply(1:length(names(out[[1]])), res_sem_time, out = out)
  colnames(out_time) <- names(out[[1]])
  t(out_time)
  
# Estaimtes Analysis ------------------------------------------------------
  
## SEM estiamtes raw data (saturated model) ##
  # Extract results per conditions
  sem_res <- lapply(seq_along(1:length(out[[1]])), # excludes problematic condtion for now
                    function(x) res_sum(out, 
                                        model = "sem", 
                                        condition = x,
                                        bias_sd = TRUE))
  
  # Show results all conditions for a given data rep
  t(sapply(1:length(sem_res),
           function(x) sem_res[[x]]$validReps))
  lapply(1:length(sem_res),
         function(x) round(sem_res[[x]]$bias_raw, 2))
  lapply(1:length(sem_res),
         function(x) sem_res[[x]]$bias_sd)
  lapply(1:length(sem_res),
         function(x) sem_res[[x]]$bias_per[-c(1:3),])
  lapply(1:length(sem_res),
         function(x) sem_res[[x]]$ci_cov)

## Linear Model: Intercept and regression coefficients ##
  lm_res <- lapply(1:length(out[[1]]),
                     function(x) res_sum(out, 
                                         model = "lm", 
                                         condition = x))
  # Show results for a given condition
  lapply(1:length(out[[1]]),
         function(x) lm_res[[x]]$bias_per)
  lapply(1:length(out[[1]]),
         function(x) lm_res[[x]]$ci_cov)

# Save Results ------------------------------------------------------------
  output <- lapply(list(sem = sem_res,
                        lm  = lm_res), 
                   function(x){
                     names(x) <- paste0("cond", seq_along(out[[1]]))
                     return(x)
                   }
  )
  
  saveRDS(
    output, 
    paste0("../output/", filename, "_res.rds")
  )
