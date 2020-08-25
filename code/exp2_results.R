### Title:    Results for experiment 2
### Porject:  Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-05-19

  rm(list = ls())
  library(xtable)
  source("./init_general.R")
  
  # Result selection
  filename <- "exp2_simOut_20200819_1743"
  
  # Read R object
  out <- readRDS(paste0("../output/", filename, ".rds"))
  
# Check presence
  out[[1]]
  names(out[[1]])
  
# Time Analyses -----------------------------------------------------------

  out_time <- sapply(1:length(names(out[[1]])), res_sem_time, out = out)
  colnames(out_time) <- names(out[[1]])
  t(out_time)
  
# Estaimtes Analysis ------------------------------------------------------
  
## SEM estiamtes raw data (saturated model) ##
  # Extract results per conditions
  semR_res <- lapply(seq_along(1:length(out[[1]]))[-2], # excludes problematic condtion for now
                     function(x) res_sum(out, 
                                         model = "semR", 
                                         condition = x))
  
  # Temporary fix for error in condition 2
  other_sem <- vector("list",8)
  other_sem[[1]] <- semR_res[[1]]
  for (i in 2:7) {
    other_sem[[i+1]] <- semR_res[[i]]
  }
  semR_res <- other_sem
  
  # Show results all conditions for a given data rep
  lapply(1:length(semR_res),
         function(x) semR_res[[x]]$bias_per)
  lapply(1:length(semR_res),
         function(x) semR_res[[x]]$ci_cov)
  
## CFA model results
  # Extract results per conditions
  CFA_res <- lapply(1:length(out[[1]]),
                     function(x) res_sum(out, 
                                         model = "CFA", 
                                         condition = x))
  
  # Results?
  lapply(1:length(out[[1]]),
         function(x) CFA_res[[x]]$bias_per[1:10, ])
  
## SEM estaimted Scored data
  # Extract results per conditions
  semS_res <- lapply(1:length(out[[1]]),
                    function(x) res_sum(out, 
                                        model = "semS", 
                                        condition = x))
  
  # Show results for a given condition
  lapply(1:length(out[[1]]),
         function(x) semS_res[[x]]$bias_per)
  lapply(1:length(out[[1]]),
         function(x) semS_res[[x]]$bias_raw)
  lapply(1:length(out[[1]]),
         function(x) semS_res[[x]]$ci_cov)

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
  output <- lapply(list(semR = semR_res,
                        CFA  = CFA_res,
                        semS = semS_res,
                        lm   = lm_res), 
                   function(x){
                     names(x) <- paste0("cond", seq_along(out[[1]]))
                     return(x)
                   }
  )

  saveRDS(
    output, 
    paste0("../output/", filename, "_res.rds") 
  )
