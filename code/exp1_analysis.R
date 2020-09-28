### Title:    Analysis of results
### Project:  Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-07-09
### Notes:    reads output form results.R script and shows the numbers that
###           are used to draw the conclusions.

  library(xtable)
  
  rm(list = ls())
  source("./init_general.R")

# Read results from a run of simulation study
  # Checks
  exp1_res <- readRDS("../output/exp1_simOut_20200731_1735_res.rds")
  sum_exp1_sem <- exp1_res$sem
  sum_exp1_lm  <- exp1_res$lm

# Bias --------------------------------------------------------------------

  # Result 1 - Easy round, all good except CART and RF
  sum_exp1_sem$cond1$cond
  sum_exp1_sem$cond1$bias_per[c(1,7,13), ]
  
  # Result 2 - Covariance trouble for Blasso and DURR w/ p 50 pm .3
  sum_exp1_sem$cond3$cond
  sum_exp1_sem$cond3$bias_per[c(1,7,13), ]
  
  # Result 3 - Variance and covariance trouble for MI-PCA and IURR, respectively
  sum_exp1_sem$cond4$cond
  sum_exp1_sem$cond4$bias_per[c(1,7,13), ]

# CI ----------------------------------------------------------------------
  
  # Result 2 - Coverage is a function of pm, not p (report mean CI coverage)
  # pm .1
  t(round(sapply(list(p50 = sum_exp1_sem$cond1$ci_cov,  # p 50 
                      p500= sum_exp1_sem$cond2$ci_cov), # p 500
               colMeans), 0))
  # p 50
  t(round(sapply(list(pm.1 = sum_exp1_sem$cond1$ci_cov,  # pm .1
                      pm.3 = sum_exp1_sem$cond3$ci_cov), # pm .3
                 colMeans), 0))
  
  # p 500
  t(round(sapply(list(pm.1 = sum_exp1_sem$cond2$ci_cov,  # pm .1
                      pm.3 = sum_exp1_sem$cond4$ci_cov), # pm .3
                 colMeans), 0))
  
# Multivairate Distance Analysis ------------------------------------------

  # Sem model estimates 
  # All parameters
  sem_ed_all <- round(t(sapply(list(p50pm.1  = sum_exp1_sem$cond1, 
                                    p500pm.1 = sum_exp1_sem$cond2, 
                                    p50pm.3  = sum_exp1_sem$cond3, 
                                    p500pm.3 = sum_exp1_sem$cond4), 
                               res_ed_est, 
                               measure = "all"
  )), 3)
  
  # Means
  sem_ed_mu <- round(t(sapply(list(p50pm.1  = sum_exp1_sem$cond1, 
                                   p500pm.1 = sum_exp1_sem$cond2, 
                                   p50pm.3  = sum_exp1_sem$cond3, 
                                   p500pm.3 = sum_exp1_sem$cond4), 
                              res_ed_est, 
                              measure = "mean"
  )), 3)
  
  # Variances
  sem_ed_var <- round(t(sapply(list(p50pm.1  = sum_exp1_sem$cond1, 
                                    p500pm.1 = sum_exp1_sem$cond2, 
                                    p50pm.3  = sum_exp1_sem$cond3, 
                                    p500pm.3 = sum_exp1_sem$cond4), 
                               res_ed_est, 
                               measure = "var"
  )), 3)
  
  # Covariances
  sem_ed_cov <- round(t(sapply(list(p50pm.1  = sum_exp1_sem$cond1, 
                                    p500pm.1 = sum_exp1_sem$cond2, 
                                    p50pm.3  = sum_exp1_sem$cond3, 
                                    p500pm.3 = sum_exp1_sem$cond4), 
                               res_ed_est, 
                               measure = "cov"
  )), 3)
  round(rbind(sem_ed_all, 
              sem_ed_mu, 
              sem_ed_var, 
              sem_ed_cov), 2)
  
  round(sem_ed_mu, 2)
  round(sem_ed_var, 2)
  round(sem_ed_cov, 2)
  
  # Make latex tables
  latex_tab_input <- rbind(sem_ed_all, rep(NA, ncol(sem_ed_all)),
                           sem_ed_mu, rep(NA, ncol(sem_ed_all)),
                           sem_ed_var, rep(NA, ncol(sem_ed_all)),
                           sem_ed_cov)
  latex_tab_input <- do.call(rbind,
                             lapply(list(sem_ed_all, 
                                         sem_ed_mu, 
                                         sem_ed_var, 
                                         sem_ed_cov), 
                                    function(x){
                                      rownames(x) <- c("p = 50, pm = .1",
                                                       "p = 50, pm = .3",
                                                       "p = 500, pm = .1",
                                                       "p = 500, pm = .3")
                                      x}
                             )
  )
  xtable(latex_tab_input,
         align = c("l", rep("c", ncol(latex_tab_input))))
  
  # or
  lapply(list(sem_ed_all, sem_ed_mu, sem_ed_var, sem_ed_cov), function(x){
    rownames(x) <- c("p = 50, pm = .1",
                     "p = 50, pm = .3",
                     "p = 500, pm = .1",
                     "p = 500, pm = .3")
    xtable(x,
           align = c("l", rep("c", ncol(x))) )
  }
  )
  

# Summary Table -----------------------------------------------------------
# This is a selected paramters version for summary paper presentation.

  col_id <- colnames(sum_exp1_sem$cond4$bias_raw)[c(12, 10, 2:9)]
  indx <- rownames(sum_exp1_sem$cond4$bias_raw)[c(1, 4, 
                                                  7, 10, 
                                                  13, 15, 25)]
  sum_exp1_sem$cond4$cond
  sum_exp1_sem$cond4$ci_cov[indx, ]
  
  sum_exp1_sem$cond4$cond
  sum_exp1_sem$cond4$bias_per[indx, ]
  
  genTableEAM <- function(x){
    # Generates section of the table for EAM turn in paper
    store <- NULL
    for (i in 1:length(indx)) {
      store <- cbind(store,
                     t(sum_exp1_sem[[x]]$bias_per[indx, col_id])[, i],
                     t(sum_exp1_sem[[x]]$ci_cov[indx, col_id])[, i])
    }
    return(store)
  }
  
  list_tables <- lapply(list(cond1=1,
                             cond2=2,
                             cond3=3,
                             cond4=4), 
                        genTableEAM)
  
  mt_table <- do.call(rbind, lapply(list_tables, rbind, NA))

  write.csv(mt_table, "~/Desktop/mt_table.csv")
  
  xtable(mt_table,
         align = c("l", rep("c", ncol(mt_table))))
  
# LM models ---------------------------------------------------------------

  # Bias
  sum_exp1_lm$cond1$cond
  sum_exp1_lm$cond1$bias_per
  sum_exp1_lm$cond3$cond
  sum_exp1_lm$cond3$bias_per
  
    # HighD
  sum_exp1_lm$cond2$cond
  sum_exp1_lm$cond2$bias_per
  sum_exp1_lm$cond4$cond
  sum_exp1_lm$cond4$bias_per
  
  # CI
  sum_exp1_lm$cond1$cond
  sum_exp1_lm$cond1$ci_cov
  sum_exp1_lm$cond3$cond
  sum_exp1_lm$cond3$ci_cov
  
    # HighD
  sum_exp1_lm$cond2$cond
  sum_exp1_lm$cond2$ci_cov
  sum_exp1_lm$cond4$cond
  sum_exp1_lm$cond4$ci_cov
  