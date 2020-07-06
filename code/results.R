### Title:    Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-05-19

rm(list = ls())

source("./init_general.R")

# out <- readRDS("../output/sim_res_20200630_1954.rds")
out <- readRDS("../output/sim_res_20200704_1612.rds")

# Time Analyses -----------------------------------------------------------

  res_sem_time(out, 1)

  out_time <- sapply(1:length(names(out[[1]])), res_sem_time, out = out)
  colnames(out_time) <- names(out[[1]])
  t(out_time)

  select_cond <- names(out[[1]])[4] # available conditions

  # Time
  res_time <- NULL
  for (i in 1:out$parms$dt_rep) {
    res_time <- rbind(res_time, out[[i]][[select_cond]]$run_time_min)
  }
  round(colMeans(res_time), 3)

# Univariate Analyses -----------------------------------------------------
  
## MLE estiamtes (saturated sem model) ##
  
  # Extract results per conditions
  out_cond1 <- res_sem_sum(out, condition = 1)
  out_cond2 <- res_sem_sum(out, condition = 2)
  out_cond3 <- res_sem_sum(out, condition = 3)
  out_cond4 <- res_sem_sum(out, condition = 4)

  # Show results for a given condition
  cnd <- 4
  names(out[[1]])[cnd]
  res_sem_sum(out, condition = cnd)$MCMC_est
  res_sem_sum(out, condition = cnd)$bias_raw
  res_sem_sum(out, condition = cnd)$bias_per
  res_sem_sum(out, condition = cnd)$ci_cov
  
## Linear Model: Intercept and regression coefficients ##
  cnd <- 4
  names(out[[1]])[cnd]
  lm_cond1 <- res_lm_sum(out, condition = 1)
  lm_cond2 <- res_lm_sum(out, condition = 2)
  lm_cond3 <- res_lm_sum(out, condition = 3)
  lm_cond4 <- res_lm_sum(out, condition = 4)
  
  lm_cond1$bias_per
  lm_cond2$bias_per
  lm_cond3$bias_per
  lm_cond4$bias_per
  
  lm_cond1$ci_cov
  lm_cond2$ci_cov
  lm_cond3$ci_cov
  lm_cond4$ci_cov
  
# Multivariate Analyses ---------------------------------------------------

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
  cov_mat <- var(rbind(avg_sts_list[[1]][,1], psd_tr_list[[1]]))
  mahalanobis(avg_sts_list[[1]][,1], psd_tr_list[[1]], cov_mat)
    # cannot be used because if methods go well, then the difference
    # between the two vectors is quite small, and the inverse of the
    # covariance matrix of highly correlated vectors is not duable
  solve(cov_mat)


# Save Results ------------------------------------------------------------
  saveRDS(
    list(cond1 = out_cond1,
         cond2 = out_cond2,
         cond3 = out_cond3,
         cond4 = out_cond4,
         parms = out$parms),
   "../output/sum_exp1.rds" 
  )
  
