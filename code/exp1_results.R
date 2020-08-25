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
                "exp1_simOut_20200801_1620")[2]
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
  # out_cond1 <- res_sem_sum(out, condition = 1)
  # out_cond2 <- res_sem_sum(out, condition = 2)
  # out_cond3 <- res_sem_sum(out, condition = 3)
  # out_cond4 <- res_sem_sum(out, condition = 4)
  # Show results for a given condition
  lapply(1:length(out[[1]]),
         function(x) sem_res[[x]]$bias_per)
  
  # cnd <- 1
  # names(out[[1]])[cnd]
  # res_sem_sum(out, condition = cnd)$MCMC_est
  # res_sem_sum(out, condition = cnd)$bias_raw
  # res_sem_sum(out, condition = cnd)$bias_per
  # res_sem_sum(out, condition = cnd)$ci_cov
  
## Linear Model: Intercept and regression coefficients ##
  # lm_cond1 <- res_lm_sum(out, condition = 1)
  # lm_cond2 <- res_lm_sum(out, condition = 2)
  # lm_cond3 <- res_lm_sum(out, condition = 3)
  # lm_cond4 <- res_lm_sum(out, condition = 4)
  
  lm_res <- lapply(1:length(out[[1]]),
                   function(x) res_sum(out, 
                                       model = "lm", 
                                       condition = x))
  
  # Show results for a given condition
  lapply(1:length(out[[1]]),
         function(x) lm_res[[x]]$bias_per)
  # cnd <- 4
  # names(out[[1]])[cnd]
  # lm_cond1$cond
  # lm_cond1$bias_per
  # lm_cond3$cond
  # lm_cond3$bias_per
  # lm_cond2$cond
  # lm_cond2$bias_per
  # lm_cond4$cond
  # lm_cond4$bias_per
  # 
  # lm_cond1$ci_cov
  # lm_cond2$ci_cov
  # lm_cond3$ci_cov
  # lm_cond4$ci_cov
  
# Multivariate Analyses ---------------------------------------------------
  
  # print(xtable(sem_ed_out_est,
  #              caption = "Experiment 1 conditions ($n = 200$)",
  #              align = c("l", rep("c", ncol(sem_ed_out_est))) ))
  # 
  # # Confidence Intervals
  # round(t(sapply(list(p50pm.1  = out_cond1, 
  #                     p50pm.3  = out_cond3, 
  #                     p500pm.1 = out_cond2, 
  #                     p500pm.3 = out_cond4), 
  #                res_ed_ci, 
  #                measure = "cov"
  # )), 1)
  # 
  # # LM model estiamtes
  # round(t(sapply(list(p50pm.1  = lm_cond1, 
  #                     p50pm.3  = lm_cond3, 
  #                     p500pm.1 = lm_cond2, 
  #                     p500pm.3 = lm_cond4), 
  #                res_ed_est, 
  #                measure = "rc"
  # )), 3)
  
# Save Results ------------------------------------------------------------
  
  output <- lapply(list(sem   = sem_res,
                        lm    = lm_res,
                        parms = out$parms), 
                   function(x){
                     names(x) <- paste0("cond", seq_along(out[[1]]))
                     return(x)
                   }
  )
  
  saveRDS(
    output, 
    paste0("../output/", filename, "_res.rds") 
  )
