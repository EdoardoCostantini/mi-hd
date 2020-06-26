### Title:    Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-05-19

rm(list = ls())

source("./functions.R")
source("./init.R")

# Current simulation results ----------------------------------------------
# If you run the simulation script you can use this directly to quick check
# Time
res_time <- NULL
for (i in 1:out[[1]]$cond_200_4$parms$dt_rep) {
  res_time <- rbind(res_time, out[[i]]$cond_200_4$run_time_min)
}
round(colMeans(res_time), 3)

# Compute Bias in a condition ---------------------------------------------

  select_cond <- names(out[[1]])[1] # available conditions

# Step 1. Obtain "true" comparison values

  full_dat_est <- NULL
  for (i in 1:out$parms$dt_rep) {
    dat <- out[[i]][[select_cond]]$dat_full
    fit <- sem(parms$lav_model, 
               data = dat,
               likelihood = "wishart")
    full_dat_est <- rbind(full_dat_est, 
                          parameterEstimates(fit)$est)
  }
  
  psd_tr_vec <- colMeans(full_dat_est) # pseudo true values

# Step 2. Compute averages of statistics

  # Row index for type of paramter
  avg_indx <- 1:out$parms$zm_n
  var_indx <- (out$parms$zm_n+1):(out$parms$zm_n*2)
  cov_indx <- (tail(var_indx, 1)+1):nrow(out[[i]][[select_cond]]$all_EST)
  
  # Store objects
  sum_stats <- matrix(0, 
                      nrow = nrow(out[[i]][[select_cond]]$all_EST), 
                      ncol = length(out$parms$methods))
  
  # Compute averages of the statistics
  for (i in 1:out$parms$dt_rep) {
      sum_stats <- sum_stats + out[[i]][[select_cond]]$all_EST
  }
  
  avg_stats <- sum_stats / out$parms$dt_rep

# Step 3. Obtain Bias
  # Raw bias
  bias <- avg_stats - psd_tr_vec
  round(bias, 3)
  
  # Bias as percentage of true value
  round(bias/psd_tr_vec*100, 3)

  # Multivariate distance
  # Euclidian distance
  
  avg_sts_list <- list(mv = avg_stats[avg_indx, ],
                       va = avg_stats[var_indx, ],
                       cv = avg_stats[cov_indx, ])
  ecd_sts_list <- t(data.frame(mv = rep(NA, length(out$parms$methods)),
                             va = rep(NA, length(out$parms$methods)),
                             cv = rep(NA, length(out$parms$methods))))
    colnames(ecd_sts_list) <- out$parms$methods
  psd_tr_list <- list(mv = psd_tr_vec[avg_indx],
                      va = psd_tr_vec[var_indx],
                      cv = psd_tr_vec[cov_indx])
  
  for (p in 1:length(avg_sts_list)) {
    for (m in 1:ncol(avg_stats)) {
      ecd_sts_list[p, m] <- dist(
        rbind(
          avg_sts_list[[p]][, m], 
          psd_tr_list[[p]]
        )
      )
    }
    }
  # Some covariances are smaller. So this is not a great measure for it
  # Or I could have all covariances here be the same
  round(ecd_sts_list, 3)

  # Mahalanobis distance 
  cov_mat <- var(rbind(bias_prct[[1]][,1], psd_tr_list[[1]]))
  mahalanobis(avg_stats[[1]][,1], psd_tr_list[[1]], cov_mat)
    # cannot be used because if methods go well, then the difference
    # between the two vectors is quite small, and the inverse of the
    # covariance matrix of highly correlated vectors is not duable
  solve(cov_mat)

# Coverage ----------------------------------------------------------------
  
  # # Row index for type of paramter
  # avg_indx <- 1:out$parms$zm_n
  # var_indx <- (out$parms$zm_n+1):(out$parms$zm_n*2)
  # cov_indx <- (tail(var_indx, 1)+1):nrow(out[[i]][[select_cond]]$all_EST)
  # psd_tr_list <- list(mv = psd_tr_vec[avg_indx],
  #                     va = psd_tr_vec[var_indx],
  #                     cv = psd_tr_vec[cov_indx])
  
  # Store objects
  str_thrs <- nrow(out[[1]][[select_cond]]$all_CI)/2 # storing threshold
  sum_stats <- matrix(0, 
                      nrow = nrow(out[[i]][[select_cond]]$all_EST), 
                      ncol = length(out$parms$methods))
  
  # Compute averages of the statistics
  for (i in 1:out$parms$dt_rep) {
    cond_est <- out[[i]][[select_cond]]$all_EST
    cond_CI <- out[[i]][[select_cond]]$all_CI
      ci_low <- cond_CI[1:str_thrs, ]
      ci_hig <- cond_CI[-(1:str_thrs), ]
    
    # General
    sum_stats <- sum_stats + sapply(1:length(out$parms$methods), 
                                    function(x){
                                      ci_out <- ci_low[, x] < psd_tr_vec &
                                        psd_tr_vec < ci_hig[, x]
                                    }
    )
  }
  
  ci_coverage <- sum_stats / out$parms$dt_rep


# Convergence -------------------------------------------------------------

  # Mice TRUE
  plot(out[[1]][[select_cond]]$imp_values$MICE_TR)
  
  # Choose imputation methods
  imp_meth <- out$parms$methods[1]
  # Plot imputation from 1st chain
  chain1 <- 1
  par(mfrow = c(3, 3))
  for (m in 1:length(imp_meth)) {
    for (v in 1:length(parms$z_m_id)) {
      imps_4plot <- out[[dt_rep]][[select_cond]]$imp_values[[imp_meth[m]]][[chain1]][[v]]
      mean_imp <- rowMeans(imps_4plot)
      plot(1:parms$iters, mean_imp, type = "l",
           main = paste0(imp_meth[m], " Mean Imputations"),
           ylim = c(mean(mean_imp)-4*sd(mean_imp), mean(mean_imp)+4*sd(mean_imp)),
           ylab = paste0("z", v), xlab = "Iteration")
  # And add imputations from other chains
      for (i in 2:(parms$chains)) {
        imps_4plot <- out[[dt_rep]][[select_cond]]$imp_values[[imp_meth[m]]][[i]][[v]]
        mean_imp <- rowMeans(imps_4plot)
        lines(1:parms$iters, mean_imp)
      }
    }
  }
