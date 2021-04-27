### Title:    Results for experiment 5
### Porject:  Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2021-04-24

  rm(list = ls())
  source("./init_general.R")
  
## Single Run  
  # Result selection
  filename <- "exp5_simOut_20210201_1850"
  
  # Read R object
  out <- readRDS(paste0("../output/", filename, ".rds"))
  
# Time Analyses -----------------------------------------------------------

  out_time <- sapply(1:length(names(out[[1]])), res_sem_time, out = out)
  colnames(out_time) <- names(out[[1]])
  t(out_time)
  
# Extract Estimates ------------------------------------------------------
  
## SEM estiamtes raw data (saturated model) ##
  # Extract results per conditions
  semR_res <- lapply(seq_along(1:length(out[[1]])),
                     function(x) res_sum(out, 
                                         model = "semR", 
                                         condition = x,
                                         bias_sd = TRUE))
  
  # Show results all conditions for a given data rep
  t(sapply(1:length(semR_res),
         function(x) semR_res[[x]]$validReps))
  lapply(1:length(semR_res),
         function(x) round(semR_res[[x]]$bias_raw[1:10,],2))
  lapply(1:length(semR_res),
         function(x) semR_res[[x]]$bias_sd)
  lapply(1:length(semR_res),
         function(x) semR_res[[x]]$bias_per[c(1:10),])
  lapply(1:length(semR_res),
         function(x) semR_res[[x]]$bias_per[-c(1:10),])
  lapply(1:length(semR_res),
         function(x) semR_res[[x]]$ci_cov)
  
## CFA model results
  # Extract results per conditions
  CFA_res <- lapply(1:length(out[[1]]),
                     function(x) res_sum(out, 
                                         model = "CFA", 
                                         condition = x))
  res_sum(out, 
          model = "CFA", 
          condition = 1)
  res_sum(out, 
          model = "CFA", 
          condition = 8)
  
  # Results?
  lapply(1:length(out[[1]]),
         function(x) CFA_res[[x]]$bias_per[1:10, ])
  
## SEM estaimted Scored data
  # Extract results per conditions
  semS_res <- lapply(1:length(out[[1]]),
                    function(x) res_sum(out, 
                                        model = "semS", 
                                        condition = x,
                                        bias_sd = TRUE))
  
  # Show results for a given condition
  lapply(1:length(out[[1]]),
         function(x) semS_res[[x]]$bias_raw[1:2,])
  lapply(1:length(out[[1]]),
         function(x) semS_res[[x]]$bias_sd)
  lapply(1:length(out[[1]]),
         function(x) semS_res[[x]]$bias_per[-c(1:2),])
  lapply(1:length(out[[1]]),
         function(x) semS_res[[x]]$ci_cov)

# Save Results ------------------------------------------------------------
  output <- lapply(list(semR = semR_res,
                        CFA  = CFA_res,
                        semS = semS_res), 
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


# Extra Details Study -----------------------------------------------------

# > Standard Deviation of estimates
  output$semR$cond4$MCMC_est[1:10, ]
  output$semR$cond4$bias_per[1:10, ]
  output$semR$cond4$bias_sd[1:10, ]
  data.frame(
    bias_per = round(colMeans(output$semR$cond4$bias_per[1:10, ]), 3)[-1],
    bias_sd = round(colMeans(output$semR$cond4$bias_sd[1:10, ]), 3),
    var_est = round(colMeans(output$semR$cond4$var_est[1:10, ]), 3)
  )
  # The average variance of the parmaters estiamtes is almost the same for the
  # different methods. In particular, bridge does not show greater variance than
  # the other methods.