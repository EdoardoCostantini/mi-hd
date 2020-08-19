### Title:    Results for experiment 2
### Porject:  Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-05-19

  rm(list = ls())
  library(xtable)
  source("./init_general.R")
  
  # Result selection
  filename <- "exp2_simOut_20200812_1449"
  
  # Read R object
  out <- readRDS(paste0("../output/", filename, ".rds"))
  
# Check presence
  out[[1]]
  
# Time Analyses -----------------------------------------------------------

  out_time <- sapply(1:length(names(out[[1]])), res_sem_time, out = out)
  colnames(out_time) <- names(out[[1]])
  t(out_time)
  
# Estaimtes Analysis ------------------------------------------------------

## SEM estiamtes raw data (saturated model) ##
  # Extract results per conditions
  semR_res <- lapply(1:length(out[[1]]),
                     function(x) res_sum(out, 
                                         model = "semR", 
                                         condition = x))
  
  # Show results for a given condition
  lapply(1:length(out[[1]]),
         function(x) print(semR_res[[x]]$bias_per))
  
## CFA model results
  # Extract results per conditions
  CFA_res <- lapply(1:length(out[[1]]),
                     function(x) res_sum(out, 
                                         model = "CFA", 
                                         condition = x))
  
  # Results?
  lapply(1:length(out[[1]]),
         function(x) print(CFA_res[[x]]$bias_per[1:10, ]))
  
## SEM estaimted Scored data
  # Extract results per conditions
  semS_res <- lapply(1:length(out[[1]]),
                    function(x) res_sum(out, 
                                        model = "semS", 
                                        condition = x))
  
  # Show results for a given condition
  lapply(1:length(out[[1]]),
         function(x) print(semS_res[[x]]$bias_per))

## Linear Model: Intercept and regression coefficients ##
  lm_res <- lapply(1:length(out[[1]]),
                     function(x) res_sum(out, 
                                         model = "lm", 
                                         condition = x))
  
  # Show results for a given condition
  lapply(1:length(out[[1]]),
         function(x) print(lm_res[[x]]$bias_per))

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
