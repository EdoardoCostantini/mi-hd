### Title:    Replication Burgette Reiter 2010 - Compare methods
### Author:   Edoardo Costantini
### Created:  2020-02-10
### Modified: 2020-02-14

  source("./functions.R")
  source("./init.R")
  
  out <- readRDS("../output/pooled_BR2010-20200214_1651.rds") # CART w/ bb
  
# Bias --------------------------------------------------------------------
  
  result_bias <- data.frame(CART = avg_bias(out, 1, parms),
                            mice = avg_bias(out, 2, parms))

# Coverage ----------------------------------------------------------------
  CART_indx <- c(1, 2)
  mi_indx <- c(3, 4)
  
  ci_cov_cart <- avg_ci_cov(out, CART_indx, parms)
  ci_cov_mice <- avg_ci_cov(out, mi_indx, parms)
  
  results_ci <- data.frame(ci_out_CART = ci_cov_cart, 
                           ci_out_mi = ci_cov_mice)
    row.names(results_ci) <- c(paste0("b", c( (1:length(parms$b_true))-1 )), "tot")
  round(results_ci,3)
  
# RMSE --------------------------------------------------------------------
  
  result_rmse <- data.frame(CART = avg_rmse(out, 1, parms = parms),
                            MICE = avg_rmse(out, 2, parms = parms))
  result_rmse
  
# Missing data descriptives
  
  missDesc <- avg_miss(out, parms = parms)

# Results -----------------------------------------------------------------

  list(Bias = result_bias,
       CI_cov = round(results_ci,3),
       rmse = result_rmse,
       miss_desc = missDesc)
  