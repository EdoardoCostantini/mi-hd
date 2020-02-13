### Title:    Replication Burgette Reiter 2010 - Compare methods
### Author:   Edoardo Costantini
### Created:  2020-JAN-10
### Modified: 2020-JAN-10
### Notes:    Compare results of replication

  output_List <- readRDS("../output/imp_out_2020-02-13_rep_1000.rds")

  source("./functions.R")                  # from source file location
  source("../../../R/functions_allpurp.R") # for bbsample function
  
# Bias --------------------------------------------------------------------
  
  B_true <- c(0, 0, .5,.5,.5,.5,.5,1,1)
  
  bias_out_list <- lapply(output_List[[1]], bias_est, x_true = B_true)
  bias_out_matr <- do.call(cbind, bias_out_list)
  cl_ind_CART <- colnames(bias_out_matr) == "pool_CART_est"
  cl_ind_mi <- colnames(bias_out_matr) == "pool_norm_est"
  
  result_bias <- data.frame(CART = round(rowMeans(bias_out_matr[,cl_ind_CART]), 3),
                            mi = round(rowMeans(bias_out_matr[,cl_ind_mi]), 3))
    row.names(result_bias) <- paste0("b", c( (1:length(B_true))-1 ) )
  result_bias

# Coverage ----------------------------------------------------------------

  CART_indx <- 1
  mi_indx <- 2
  
  ci_out_mi <- ci_out_CART <- 0

  # For CART
  for (dt in 1:length(output_List[[1]])) {
    confi <- output_List[[3]][[dt]][[CART_indx]]
    ci_out_CART <- ci_out_CART + as.numeric(c(B_true > confi[,1] & B_true < confi[,2]))
  }
  
  bw_cov_cart <- ci_out_CART/length(output_List[[1]])
  tot_cov_cart <- sum(ci_out_CART)/(length(B_true)*length(output_List[[1]]))
  
  # For MI
  for (dt in 1:length(output_List[[1]])) {
    confi <- output_List[[3]][[dt]][[mi_indx]]
    ci_out_mi <- ci_out_mi + as.numeric(c(B_true > confi[,1] & B_true < confi[,2]))
  }
  
  bw_cov_mice <- ci_out_mi/length(output_List[[1]])
  tot_cov_mice <- sum(ci_out_mi)/(length(B_true)*length(output_List[[1]]))
  
  results_ci <- data.frame(ci_out_CART = c(bw_cov_cart, tot_cov_cart), 
                           ci_out_mi = c(bw_cov_mice, tot_cov_mice))
    row.names(results_ci) <- c(paste0("b", c( (1:length(B_true))-1 )), "tot")
  round(results_ci,3)
  
# RMSE --------------------------------------------------------------------
  
  list_rmse <- lapply(output_List[[2]], rmse_est)

  rmse_out_matr <- do.call(cbind, list_rmse)
  cl_ind_CART <- colnames(rmse_out_matr) == "pool_CART_est"
  cl_ind_mi <- colnames(rmse_out_matr) == "pool_mi_est"
  
  result_rmse <- round(rowMeans(rmse_out_matr), 3)[-3]
  result_rmse

# Results -----------------------------------------------------------------

  list(Bias = result_bias,
       CI_cov = round(results_ci,3),
       rmse = result_rmse)
  