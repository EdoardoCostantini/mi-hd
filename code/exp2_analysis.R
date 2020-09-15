### Title:    Analysis of results from experiment 2
### Project:  Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-07-09
### Notes:    reads output form results.R script and shows the numbers that
###           are used to draw the conclusions.

  library(xtable)
  
  rm(list = ls())
  source("./init_general.R")

# Read results from a run of simulation study
  # exp2_res <- readRDS("../output/exp2_simOut_20200812_1449_res.rds") # way out
  filename <- "exp2_simOut_20200819_1743"
  out <- readRDS(paste0("../output/", filename, ".rds"))
  exp2_res <- readRDS(paste0("../output/", filename, "_res.rds"))

# Bias --------------------------------------------------------------------
  
# Recap of set up  

  out$conds
  
  # c  lv  pm   fl ridge
  # 1  10 0.1 high 1e-01
  # 2 100 0.1 high 1e-07
  # 3  10 0.3 high 1e-01
  # 4 100 0.3 high 1e-07
  # ------------------ #
  # 5  10 0.1  low 1e-01
  # 6 100 0.1  low 1e-07
  # 7  10 0.3  low 1e-01
  # 8 100 0.3  low 1e-07
  
  list(
    it_number = out$parms$n_it,
    mis_var   = out$parms$z_m_id,
    miss_type = out$parms$missType,
    lv_rm_x   = out$parms$rm_x,
    lv_number = "condition specific: 10 or 100",
    factor_loadings = "a) runif btw .9 and .97; b) runif btw .5 and .6",
    substantive = "SEM, CFA raw data; SEM, LM scored data (mean of items)"
  )
  
  # Condition indexes
  cindex_lh <- c(1:4)
  cindex_hd <- c(2, 4)
  cindex_hp <- c(3, 4)
  
#> SEM scored ####

  # Var Covar
  # PCA issue with variances remains but good performances
  # HIGH FACTOR LOADINGS
  lapply(exp2_res$semS,
         function(x) x$bias_per[-c(1:2),])[cindex_lh]
  # LOW FACTOR LOADINGS
  lapply(exp2_res$semS,
         function(x) x$bias_per[-c(1:2),])[-cindex_lh]
  
  # MEANS in terms of standard deviations from the reference
  lapply(exp2_res$semS,
         function(x) x$bias_sd)[c(5:8)]
  
  # CONFIDENCE INTERVALS
  # Look at sizes
  lapply(exp2_res$semS,
         function(x) x$ci_cov)[c(4,8)]
  # Look at ED
  t(sapply(exp2_res$semS,
           res_ed_ci))
  
#> SEM raw data ####
  # t(sapply(exp2_res$semR,
  #          function(x) x$validReps))
  
  # EUCLIDEAN distance measure for means
  # LOW FACTOR LOADINGS
  t(sapply(exp2_res$semR,
           res_ed_est, index = 1:10))[5:8, ]
  
  # VARIANCES
  # Good PCA good for all but highest condition with low factor loadings
  # Last condition is the interesting one: PCA and IURR perform well in all
  # other conditions but when factor loadings are smaller then PCA starts to
  # show its biased variances problem, and IURR its biased covariances problem
  lapply(exp2_res$semR,
         function(x) x$bias_per[-c(1:10),])[c(4,8)]
  
  # PCA variance struggle
  t(sapply(exp2_res$semR,
           res_ed_est, index = 11:20))
  
  # PCA covariances domination
  t(sapply(exp2_res$semR,
           res_ed_est, index = -c(1:20)))
  
  # Confidence Intervals one by one too crowded
  # lapply(exp2_res$semR,
  #        function(x) x$ci_cov)[c(4, 8)]
  
  # Look at Euclidean distance measure
  # CI are great for PCA!
  t(sapply(exp2_res$semR,
           res_ed_ci))[cindex_lh, ]
  t(sapply(exp2_res$semR,
           res_ed_ci))[-cindex_lh, ]
  
#> CFA raw data ####
  # FACTOR LOADINGS: PCA does not even budge with high dim and low factor 
  # loadings
  lapply(exp2_res$CFA,
         function(x) x$bias_per[1:10, ])[cindex_lh]
  lapply(exp2_res$CFA,
         function(x) x$bias_per[1:10, ])[-cindex_lh]
  # # CONFIDENCE INTERVALS ARE GOOD AS WELL
  # lapply(exp2_res$CFA,
  #        function(x) x$ci_cov[1:10, ])[cindex_lh]
  # lapply(exp2_res$CFA,
  #        function(x) x$ci_cov[1:10, ])[-cindex_lh]
  

  
#> lm scored data ####
  
  # # Bias
  # # High Factor loadings
  # lapply(exp2_res$lm,
  #        function(x) x$bias_per)[cindex_lh]
  # # Low factor loadings
  # lapply(exp2_res$lm,
  #        function(x) x$bias_per)[-cindex_lh]
  # 
  # # CI
  # lapply(exp2_res$lm,
  #        function(x) x$ci_cov)[cindex_lh]
  # lapply(exp2_res$lm,
  #        function(x) x$ci_cov)[-cindex_lh]

# # Summary Table -----------------------------------------------------------
# # This is a selected paramters version for summary paper presentation.
# 
#   col_id <- colnames(sum_exp1_sem$cond4$bias_raw)[c(12, 10, 2:9)]
#   indx   <- rownames(sum_exp1_sem$cond4$bias_raw)[c(1, 4, 
#                                                     7, 10, 
#                                                     13, 15, 25)]
#   sum_exp1_sem$cond4$cond
#   sum_exp1_sem$cond4$ci_cov[indx, ]
#   
#   sum_exp1_sem$cond4$cond
#   sum_exp1_sem$cond4$bias_per[indx, ]
#   
#   genTableEAM <- function(x){
#     # Generates section of the table for EAM turn in paper
#     store <- NULL
#     for (i in 1:length(indx)) {
#       store <- cbind(store,
#                      t(sum_exp1_sem[[x]]$bias_per[indx, col_id])[, i],
#                      t(sum_exp1_sem[[x]]$ci_cov[indx, col_id])[, i])
#     }
#     return(store)
#   }
#   
#   list_tables <- lapply(list(cond1=1,
#                              cond2=2,
#                              cond3=3,
#                              cond4=4), 
#                         genTableEAM)
#   
#   mt_table <- do.call(rbind, lapply(list_tables, rbind, NA))
# 
#   write.csv(mt_table, paste0("../output/", filename, "_table.csv") )
