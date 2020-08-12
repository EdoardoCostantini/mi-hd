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
  
# Time Analyses -----------------------------------------------------------

  out_time <- sapply(1:length(names(out[[1]])), res_sem_time, out = out)
  colnames(out_time) <- names(out[[1]])
  t(out_time)
  
# Univariate Analyses -----------------------------------------------------
  names(out[[1]]$`cond_10_0.1_high_1e-05`)
  c("semR","CFA","semS")
  
## SEM estiamtes raw data (saturated model) ##
  # Extract results per conditions
  semR_cond1 <- res_sum(out, model = "semR", condition = 1)
  semR_cond2 <- res_sum(out, model = "semR", condition = 2)
  semR_cond3 <- res_sum(out, model = "semR", condition = 3)
  semR_cond4 <- res_sum(out, model = "semR", condition = 4)
  
  # Show results for a given condition
  semR_cond1$MCMC_est
  semR_cond1$bias_raw
  semR_cond1$bias_per
  semR_cond1$ci_cov
  
## CFA model results
  # Extract results per conditions
  CFA_cond1 <- res_sum(out, model = "CFA", condition = 1)
  CFA_cond2 <- res_sum(out, model = "CFA", condition = 2)
  CFA_cond3 <- res_sum(out, model = "CFA", condition = 3)
  CFA_cond4 <- res_sum(out, model = "CFA", condition = 4)
  
  # Results?
  CFA_cond1$bias_per[1:10, ]
  CFA_cond2$bias_per[1:10, ]
  CFA_cond3$bias_per[1:10, ]
  CFA_cond4$bias_per[1:10, ]
  
## SEM estaimted Scored data
  # Extract results per conditions
  semS_cond1 <- res_sum(out, model = "semS", condition = 1)
  semS_cond2 <- res_sum(out, model = "semS", condition = 2)
  semS_cond3 <- res_sum(out, model = "semS", condition = 3)
  semS_cond4 <- res_sum(out, model = "semS", condition = 4)
  
  # Show results for a given condition
  semS_cond1$bias_per
  semS_cond2$bias_per
  semS_cond3$bias_per
  semS_cond4$bias_per

## Linear Model: Intercept and regression coefficients ##
  lm_cond1 <- res_sum(out, model = "lm", condition = 1)
  lm_cond2 <- res_sum(out, model = "lm", condition = 2)
  lm_cond3 <- res_sum(out, model = "lm", condition = 3)
  lm_cond4 <- res_sum(out, model = "lm", condition = 4)
  
  # Show results for a given condition
  lm_cond1$bias_per
  lm_cond3$bias_per
  lm_cond2$bias_per
  lm_cond4$bias_per

# Save Results ------------------------------------------------------------
  saveRDS(
    list(semR = list(cond1 = semR_cond1,
                     cond2 = semR_cond2,
                     cond3 = semR_cond3,
                     cond4 = semR_cond4),
         CFA = list(cond1 = CFA_cond1,
                    cond2 = CFA_cond2,
                    cond3 = CFA_cond3,
                    cond4 = CFA_cond4),
         semS = list(cond1 = semS_cond1,
                     cond2 = semS_cond2,
                     cond3 = semS_cond3,
                     cond4 = semS_cond4),
         lm = list(cond1 = lm_cond1,
                   cond2 = lm_cond2,
                   cond3 = lm_cond3,
                   cond4 = lm_cond4)), 
    paste0("../output/", filename, "_res.rds") 
  )
