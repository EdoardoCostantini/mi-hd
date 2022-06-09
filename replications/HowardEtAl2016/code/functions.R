### Title:    Replication Howard et al 2015
### Author:   Anonymized for peer review
### Created:  2020-05-18

# data generation ---------------------------------------------------------

genData <- function(cnd, parms){
  # For internals
  # cnd <- conds[1,]
  
  Sigma <- matrix(parms$cor_fix, 
                  nrow = parms$n_vars, ncol = parms$n_vars)
    colnames(Sigma) <- rownames(Sigma) <- parms$vars
  diag(Sigma) <- 1
  Sigma["Y", paste0("A", seq(1, 8))] <- cnd["YA_rho"]
  Sigma[paste0("A", seq(1, 8)), "Y"] <- cnd["YA_rho"]
  
  # Create population dataset
  dt_pop <- as.data.frame(rmvnorm(parms$n, 
                                  mean = rep(parms$mean, parms$n_vars),
                                  sigma = Sigma))
  colnames(dt_pop) <- parms$vars
  
  # Sample from population
  idx <- base::sample(1:parms$n, cnd["N"], replace = FALSE)
  dt_smp <- dt_pop[idx, ]
  
  # Return data  
  return( list(pop = dt_pop,
               smp = dt_smp) )
}

# Example use
# set.seed(1234)
# dt <- genData(conds[1, ], parms)
# colMeans(dt$smp)
# cov(dt$smp)
# cor(dt$smp)

# lm(y ~ z1 + z2 + z3 + z4 + z5, data = Xy) # true y model

impose_NL_Miss <- function(dt_full, cnd, parms){
  # Given a fully observed dataset and param object containing the regression
  # coefficients of the model to impose missingness, it returns a version of
  # the original data with imposed missingness on z1, and its missingness
  # idicator (R)
  
  # cnd <- conds[1,]
  # dt_full <- dt$smp
  
  A1A2 <- dt_full$A1*dt_full$A2 # interaction score
  mis_idx <- A1A2 <= sort(A1A2)[floor(cnd["N"] * cnd["miss_r"])]
  dt_full[mis_idx, "Y"] <- NA
  
  return(dt_full)
}

impose_L_Miss <- function(dt_full, cnd, parms){
  # Given a fully observed dataset and param object containing the regression
  # coefficients of the model to impose missingness, it returns a version of
  # the original data with imposed missingness on z1, and its missingness
  # idicator (R)
  
  # cnd <- conds[1,]
  # dt_full <- dt$smp
  
  mis_idx <- dt_full$A1 <= sort(dt_full$A1)[floor(cnd["N"] * cnd["miss_r"])]
  dt_full[mis_idx, "Y"] <- NA
  
  return(dt_full)
}

impose_MCAR <- function(dt_full, cnd, parms){
  # cnd <- conds[1, ]
  for(i in 1:8){
    mis_idx <- sample(1:cnd["N"], floor(cnd["N"]*cnd["miss_r"]))
    dt_full[mis_idx, parms$vars[-c(1,2)][i] ] <- NA
  }
  return(dt_full)
}

# Results -----------------------------------------------------------------
bias <- function(sem_par_out, parms){
  
  # sem_par_out <- pe_GSp
  # Select parameters of itnerest
  sem_par_names <- sapply(1:nrow(sem_par_out), 
                          function(x) paste0(sem_par_out[x, 1:3], collapse = ''))
  sem_par_idx <- sapply(1:length(parms$out_ref),
                        function(x) which(grepl(parms$out_ref[x], sem_par_names)))
  
  PE <- sem_par_out[sem_par_idx, 4]
  
  # Compute Bias
  Bias_res <- PE - c(parms$cor_fix, parms$mean, parms$sd**2)
  
  return(Bias_res)
}

cove <- function(sem_par_out, parms){
  
  # sem_par_out <- pe_GSp
  # Select parameters of itnerest
  sem_par_names <- sapply(1:nrow(sem_par_out), 
                          function(x) paste0(sem_par_out[x, 1:3], collapse = ''))
  sem_par_idx <- sapply(1:length(parms$out_ref),
                        function(x) which(grepl(parms$out_ref[x], sem_par_names)))
  
  PE <- sem_par_out[sem_par_idx, 4]
  CI <- sem_par_out[sem_par_idx, c(5, 6)]
  
  # Compute confidence coverage
  CI_res <- NULL
  for (i in 1:3) {
    CI_res[i] <- PE[i] > CI[i, 1] && PE[i] < CI[i, 2]
  }
  return(CI_res)
}

extract_results <- function(cond_name, output, dt_rep){
  # Example input
  # cond_name <- names(out[[1]])[2]
  # output <- out
  # dt_rep = out[[1]][[cond_name]]$parms$dt_rep
  
  # Bias
  store_sum <- vector("list", dt_rep)
  
  for (i in 1:dt_rep) {
    store_sum[[i]] <- output[[i]][[cond_name]]$rep_bias
  }
  
  bias_out <- round(Reduce("+", store_sum)/dt_rep, 3)
  rownames(bias_out) <- parms$out_ref
  
  # Average Coverage
  store_sum <- vector("list", dt_rep)
  
  for (i in 1:dt_rep) {
    store_sum[[i]] <- output[[i]][[cond_name]]$rep_cov
  }
  
  CI_out <- Reduce("+", store_sum)/dt_rep
    rownames(CI_out) <- parms$out_ref
  
  resu <- list(bias = bias_out, 
               CI = round(CI_out, 3))
  return(resu)
}
