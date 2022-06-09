### Title:    Replication Burgette Reiter 2010 - Compare methods
### Author:   Anonymized for peer review
### Created:  2020-02-10
### Modified: 2020-02-14

  source("./functions.R")
  source("./init.R")
  
  out <- readRDS("../output/pooled_BR2010-mc-20200218_2232.rds") # CART w/ bb
  
# Bias --------------------------------------------------------------------
  
  result_bias <- data.frame(CARTbb = avg_bias(out, 1, parms),
                            CARTd = avg_bias(out, 3, parms),
                            mice = avg_bias(out, 2, parms))
    rownames(result_bias) <- paste0("b", c( (0:(length(parms$b_true)-1)) ))

# Coverage ----------------------------------------------------------------
  CARTbb_indx <- c(1, 2)
  CARTd_indx <- c(5, 6)
  mi_indx <- c(3, 4)
  
  results_ci <- data.frame(CARTbb = avg_ci_cov(out, CARTbb_indx, parms), 
                           CARTd = avg_ci_cov(out, CARTd_indx, parms), 
                           mice = avg_ci_cov(out, mi_indx, parms))
  
    row.names(results_ci) <- c(paste0("b", c( (1:length(parms$b_true))-1 )), "tot")
  round(results_ci, 3)
  
# RMSE --------------------------------------------------------------------
  
  result_rmse <- data.frame(CARTbb = avg_rmse(out, 1, parms = parms),
                            CARTd = avg_rmse(out, 3, parms = parms),
                            mice = avg_rmse(out, 2, parms = parms))
  result_rmse
  
# Missing data descriptives
  
  missDesc <- avg_miss(out, parms = parms)

# Results -----------------------------------------------------------------

  list(Bias = result_bias,
       CI_cov = round(results_ci,3),
       rmse = result_rmse,
       miss_desc = missDesc)
  