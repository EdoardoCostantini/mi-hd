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
  # filename <- "sim_res_20200731_1735"
  
  # Read R object
  out <- readRDS(paste0("../output/", filename, ".rds"))

  out$parms$dt_rep
  out$parms$iters
  
# Time Analyses -----------------------------------------------------------

  out_time <- sapply(1:length(names(out[[1]])), res_sem_time, out = out)
  colnames(out_time) <- names(out[[1]])
  t(out_time)
  
# Univariate Analyses -----------------------------------------------------
  
## MLE estiamtes (saturated sem model) ##
  
  # Extract results per conditions
  out_cond1 <- res_sem_sum(out, condition = 1)
  out_cond2 <- res_sem_sum(out, condition = 2)
  out_cond3 <- res_sem_sum(out, condition = 3)
  out_cond4 <- res_sem_sum(out, condition = 4)
  
  # Show results for a given condition
  cnd <- 1
  names(out[[1]])[cnd]
  res_sem_sum(out, condition = cnd)$MCMC_est
  res_sem_sum(out, condition = cnd)$bias_raw
  res_sem_sum(out, condition = cnd)$bias_per
  res_sem_sum(out, condition = cnd)$ci_cov
  
## Linear Model: Intercept and regression coefficients ##
  lm_cond1 <- res_lm_sum(out, condition = 1)
  lm_cond2 <- res_lm_sum(out, condition = 2)
  lm_cond3 <- res_lm_sum(out, condition = 3)
  lm_cond4 <- res_lm_sum(out, condition = 4)
  
  # Show results for a given condition
  cnd <- 4
  names(out[[1]])[cnd]
  lm_cond1$cond
  lm_cond1$bias_per
  lm_cond3$cond
  lm_cond3$bias_per
  lm_cond2$cond
  lm_cond2$bias_per
  lm_cond4$cond
  lm_cond4$bias_per
  
  lm_cond1$ci_cov
  lm_cond2$ci_cov
  lm_cond3$ci_cov
  lm_cond4$ci_cov
  
# Multivariate Analyses ---------------------------------------------------
  
  print(xtable(sem_ed_out_est,
               caption = "Experiment 1 conditions ($n = 200$)",
               align = c("l", rep("c", ncol(sem_ed_out_est))) ))
  
  # Confidence Intervals
  round(t(sapply(list(p50pm.1  = out_cond1, 
                      p50pm.3  = out_cond3, 
                      p500pm.1 = out_cond2, 
                      p500pm.3 = out_cond4), 
                 res_ed_ci, 
                 measure = "cov"
  )), 1)
  
  # LM model estiamtes
  round(t(sapply(list(p50pm.1  = lm_cond1, 
                      p50pm.3  = lm_cond3, 
                      p500pm.1 = lm_cond2, 
                      p500pm.3 = lm_cond4), 
                 res_ed_est, 
                 measure = "rc"
  )), 3)
  
# Save Results ------------------------------------------------------------
  saveRDS(
    list(cond1 = out_cond1,
         cond2 = out_cond2,
         cond3 = out_cond3,
         cond4 = out_cond4,
         parms = out$parms),
    paste0("../output/", filename, "sum_exp1_sem.rds") 
  )
  
  saveRDS(
    list(cond1 = lm_cond1,
         cond2 = lm_cond2,
         cond3 = lm_cond3,
         cond4 = lm_cond4,
         parms = out$parms),
    paste0("../output/", filename, "sum_exp1_lm.rds") 
  )