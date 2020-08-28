### Title:    Analysis of results from experiment 2
### Project:  Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-08-27
### Notes:    reads output form results.R script and shows the numbers that
###           are used to draw the conclusions.

  library(xtable)
  
  rm(list = ls())
  source("./init_general.R")

# Read results from a run of simulation study
  # exp3_res <- readRDS("../output/exp2_simOut_20200812_1449_res.rds") # way out
  filename <- "exp3_simOut_20200825_0951"
  out      <- readRDS(paste0("../output/", filename, ".rds"))
  exp3_res <- readRDS(paste0("../output/", filename, "_res.rds"))


# Summary of set up -------------------------------------------------------
  out$conds
  
#     pm ridge   r2 int_sub int_rm int_da   p
# 1  0.3 1e-08 0.65   FALSE  FALSE  FALSE  50
# 2  0.3 1e-01 0.65    TRUE  FALSE  FALSE  50
# 3  0.3 1e-04 0.65   FALSE   TRUE  FALSE  50
# 4  0.3 1e-01 0.65    TRUE   TRUE  FALSE  50
# ----------------------------------------- #
# 5  0.3 1e-08 0.65   FALSE  FALSE   TRUE  50
# 6  0.3 1e-08 0.65    TRUE  FALSE   TRUE  50
# 7  0.3 1e-07 0.65   FALSE   TRUE   TRUE  50
# 8  0.3 1e-07 0.65    TRUE   TRUE   TRUE  50
# ----------------------------------------- #
# 9  0.3 1e-07 0.65   FALSE  FALSE  FALSE 500
# 10 0.3 1e-07 0.65    TRUE  FALSE  FALSE 500
# 11 0.3 1e-07 0.65   FALSE   TRUE  FALSE 500
# 12 0.3 1e-07 0.65    TRUE   TRUE  FALSE 500
# ----------------------------------------- #
  
  # Missing variables + interaction term in JAV style
  out$parms$z_m_id
  # Substantive models
  list(formula = out$parms$frm,
       formula_int = out$parms$frm_int,
       coefs = c(out$parms$b_main, out$parms$b_int))
  # Response models
  list(regular = out$parms$rm_x,
       wInter = out$parms$rm_x_int)
  
#> SEM raw data ####
  t(sapply(exp3_res$sem,
           function(x) x$validReps))
  
  # Means
  lapply(exp3_res$sem,
         function(x) x$bias_raw[1:3,]
  )
  
  # Euclidean distance measure
  exp3_res$sem[[3]]$MCMC_est[1:10, ] # rows of interest
  t(sapply(exp3_res$sem,
           res_ed_est, index = 1:10))[cindex_lh, ]
  t(sapply(exp3_res$sem,
           res_ed_est, index = 1:10))[-cindex_lh, ]
  
  # Variances
  
  # No DA
  lapply(exp3_res$sem,
         function(x) x$bias_per[-c(1:3),])[cindex_lnda]
  
  # Yes DA
  lapply(exp3_res$sem,
         function(x) x$bias_per[-c(1:3),])[cindex_lyda]
    # Last condition is the interesting one: PCA and IURR perform well in all
    # other conditions but when factor loadings are smaller then PCA starts to
    # show its biased variances problem, and IURR its biased covariances problem
  
  # Int in substantive model: no DA vs DA
  lapply(exp3_res$sem,
         function(x) x$bias_per[-c(1:3),])[c(2, 6)]
  
  # Int in response model: no DA vs DA
  lapply(exp3_res$sem,
         function(x) x$bias_per[-c(1:3),])[c(3, 7)]
  
  # Confidence Intervals
  lapply(exp3_res$sem,
         function(x) x$ci_cov)[cindex_lh]
  
  t(sapply(exp3_res$sem,
           res_ed_ci))[cindex_lh, ]
  t(sapply(exp3_res$sem,
           res_ed_ci))[-cindex_lh, ]
    # CI are great for PCA
  
#> lm scored data ####
  
  # INTERACTION PRESENCE
  # - CART and DURR on top?
  # - presence in RM does not make things harder/worst
  
  # Low dim
  lapply(exp3_res$lm,
         function(x) x$bias_per)[1:4]
  
  # High dim
  lapply(exp3_res$lm,
         function(x) x$bias_per)[9:12]
  
  # NO INTERACTIONS PCA does not well with DA
  lapply(exp3_res$lm,
         function(x) x$bias_per)[c(5:8)] # could be a PCA selection issue?
  
  # but does not suffer much from increased dimensionality
  lapply(exp3_res$lm,
         function(x) x$bias_per)[c(1,3, 9,11)]
  
  # INTERACTION MANAGMENT
  # Int in substantive model: no DA vs DA
  lapply(exp3_res$lm,
         function(x) x$bias_per)[c(2, 6)]
  
  # Int in rm model: no DA vs DA
  lapply(exp3_res$lm,
         function(x) x$bias_per)[c(3, 7)]
  
  # Int everywhere: no DA vs DA
  lapply(exp3_res$lm,
         function(x) x$bias_per)[c(4, 8)]
  
  # # Confidence inTerval patterns are similar
  # lapply(exp3_res$lm,
  #        function(x) x$ci_cov)
  
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
