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
  filename <- "exp2_simOut_20200819_1743_res"
  exp2_res <- readRDS(paste0("../output/", filename, ".rds"))

# Bias --------------------------------------------------------------------
  exp2_res$conds
  
  # Condition indexes
  cindex_lh <- c(1:4)
    runif(n_it_tot, .9, .97) #lh
    runif(n_it_tot, .5, .6)  #ll
  cindex_hd <- c(2, 4)
  cindex_hp <- c(3, 4)
  
#> SEM raw data ####
  t(sapply(exp2_res$semR,
           function(x) x$validReps))
  exp2_res$semR[[4]]$MCMC_est
  
  # Means
  lapply(exp2_res$semR,
         function(x) x$bias_raw[1:10,]
  )[cindex_lh]
  
  lapply(exp2_res$semR,
         function(x) x$bias_raw[1:10,]
         )[cindex_lh]
    # Bias in percent makes no sense for means = 0
  
  lapply(exp2_res$semR,
         function(x) x$bias_sd)[cindex_lh]
    # High factor loadings: PCA and IURR are impressivly good
  
  lapply(exp2_res$semR,
         function(x) x$bias_sd)[-cindex_lh]
    # Low factor loadings: same pattern
  
  # Euclidean distance measure
  exp2_res$semR[[3]]$MCMC_est[1:10, ] # rows of interest
  t(sapply(exp2_res$semR,
           res_ed_est, index = 1:10))[cindex_lh, ]
  t(sapply(exp2_res$semR,
           res_ed_est, index = 1:10))[-cindex_lh, ]
  
  # Variances
  lapply(exp2_res$semR,
         function(x) x$bias_per[-c(1:10),])[cindex_lh]
  lapply(exp2_res$semR,
         function(x) x$bias_per[-c(1:10),])[-cindex_lh]
    # Last condition is the interesting one: PCA and IURR perform well in all
    # other conditions but when factor loadings are smaller then PCA starts to
    # show its biased variances problem, and IURR its biased covariances problem
  
  # Confidence Intervals
  lapply(exp2_res$semR,
         function(x) x$ci_cov)[cindex_lh]
  
  t(sapply(exp2_res$semR,
           res_ed_ci))[cindex_lh, ]
  t(sapply(exp2_res$semR,
           res_ed_ci))[-cindex_lh, ]
    # CI are great for PCA
  
#> CFA raw data ####
  lapply(exp2_res$CFA,
         function(x) x$bias_per[1:10, ])[cindex_lh]
  lapply(exp2_res$CFA,
         function(x) x$bias_per[1:10, ])[-cindex_lh]
  
  lapply(exp2_res$CFA,
         function(x) x$ci_cov[1:10, ])[cindex_lh]
  lapply(exp2_res$CFA,
         function(x) x$ci_cov[1:10, ])[-cindex_lh]
  
#> SEM scored ####
  t(sapply(exp2_res$semS,
           function(x) x$validReps))
  
  # Means of meaned scored data
  lapply(exp2_res$semS,
         function(x) round(x$bias_raw[1:2,],2))[cindex_lh]
  
  lapply(exp2_res$semS,
         function(x) x$bias_sd)[cindex_lh]
  lapply(exp2_res$semS,
         function(x) x$bias_sd)[-cindex_lh]
  
  # Var Covar
  lapply(exp2_res$semS,
         function(x) x$bias_per[-c(1:2),])[cindex_lh]
  lapply(exp2_res$semS,
         function(x) x$bias_per[-c(1:2),])[-cindex_lh]
  
  # Confidence Interval
  lapply(exp2_res$semS,
         function(x) x$ci_cov)[cindex_lh]
  
  t(sapply(exp2_res$semS,
           res_ed_ci))[cindex_lh, ]
  t(sapply(exp2_res$semS,
           res_ed_ci))[-cindex_lh, ]
  
#> lm scored data ####
  
  # Bias
  lapply(exp2_res$lm,
         function(x) x$bias_per)[cindex_lh]
  lapply(exp2_res$lm,
         function(x) x$bias_per)[-cindex_lh]
  
  # CI
  lapply(exp2_res$lm,
         function(x) x$ci_cov)[cindex_lh]
  lapply(exp2_res$lm,
         function(x) x$ci_cov)[-cindex_lh]

# Summary Table -----------------------------------------------------------
# This is a selected paramters version for summary paper presentation.

  col_id <- colnames(sum_exp1_sem$cond4$bias_raw)[c(12, 10, 2:9)]
  indx   <- rownames(sum_exp1_sem$cond4$bias_raw)[c(1, 4, 
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

  write.csv(mt_table, paste0("../output/", filename, "_table.csv") )
