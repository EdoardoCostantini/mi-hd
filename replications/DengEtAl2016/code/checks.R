### Title:    Replication Zhao Long 2016 - Check certain key elements
### Author:   Anonymized for peer review
### Created:  2020-02-20
### Modified: 2020-02-20
### Notes:    Check whether some setups of the simulation are achieved

source("./init.R")
source("./functions.R")

# Missingness -------------------------------------------------------------

  cond <- conds[1,]
  obj <- NULL
  for (i in 1:parms$dt_rep) {
    Xy <- genData(cond, parms)
    Xy_mis <- imposeMiss(Xy, parms)$Xy_miss
    obj[i] <- mean(rowSums(is.na(Xy_mis)) != 0)
  }
  round(mean(obj), 1) # approximately .4
  # Should be .4 missingness on this variable on average
  

# Check no impoutation results --------------------------------------------
# Try obtaining the same results as in table 1 for the analysis without imputation
# Or simple mice imputation.
  set.seed(1234)
  
  store_coef_gs <- vector("list", 10)
  store_bias_gs <- vector("list", 10)
  store_se_gs <- vector("list", 10)
  store_coef_cc <- vector("list", 10)
  store_bias_cc <- vector("list", 10)
  store_se_cc <- vector("list", 10)
  
  for (rr in 1:10) {
    reps <- 200
    coef_GS <- matrix(rep(NA, 6*reps), ncol = 6)
    bias_GS <- matrix(rep(NA, 6*reps), ncol = 6)
    se_GS <- matrix(rep(NA, 6*reps), ncol = 6)
    CI_GS <- matrix(rep(NA, 6*reps), ncol = 6)
    
    coef_CC <- matrix(rep(NA, 6*reps), ncol = 6)
    bias_CC <- matrix(rep(NA, 6*reps), ncol = 6)
    se_CC <- matrix(rep(NA, 6*reps), ncol = 6)
    CI_CC <- matrix(rep(NA, 6*reps), ncol = 6)
    
    bias_TR <- matrix(rep(NA, 6*reps), ncol = 6)
    CI_TR <- matrix(rep(NA, 6*reps), ncol = 6)
    
    for (i in 1:reps) {
      # fully observed data analysis
      Xy <- genData(conds[1, ], parms)
      mod_GS <- lm(parms$formula, data = Xy)
      summa_GS <- summary(mod_GS)
      
      coef_GS[i, ] <- mod_GS$coefficients
      bias_GS[i, ] <- mod_GS$coefficients - 1
      se_GS[i, ] <- summa_GS$coefficients[, 2]
      CI_GS[i, ] <- check_cover(confint(mod_GS))
      
      # Complete case analysis
      Xy_mis <- imposeMiss(Xy, parms)$Xy_miss
      mod_CC <- lm(parms$formula, data = Xy_mis)
      summa_CC <- summary(mod_CC)
      
      coef_CC[i, ] <- mod_CC$coefficients
      bias_CC[i, ] <- mod_CC$coefficients - 1
      se_CC[i, ] <- summa_CC$coefficients[, 2]
      CI_CC[i, ] <- check_cover(confint(mod_CC))
    }
    store_coef_gs[[rr]] <- coef_GS
    store_bias_gs[[rr]] <- bias_GS
    store_se_gs[[rr]]   <- se_GS
    store_bias_cc[[rr]] <- bias_CC
    store_se_cc[[rr]]   <- se_CC
    store_coef_cc[[rr]] <- coef_CC
  }
  
  # Range of Mean bias
  list(
    GS = round(t(sapply(store_bias_gs, colMeans)), 3),
    CC = round(t(sapply(store_bias_cc, colMeans)), 3)
  )
  
  # Mean Standard errors
  list(
    GS = round(t(sapply(store_se_gs, colMeans)), 3),
    CC = round(t(sapply(store_se_cc, colMeans)), 3)
  )
  
  # Monte Carlo standard deviation
  list(
    GS = round(t(sapply(1:10, function(x) apply(store_coef_cc[[x]], 2, sd) )), 3),
    CC = round(t(sapply(1:10, function(x) apply(store_coef_cc[[x]], 2, sd) )), 3)
  )

  reps <- 500
  
  # LM related
  b_fl <- b_ms <- 
    matrix(NA, nrow = reps, ncol = 6)
  R2 <- matrix(NA, nrow = reps, ncol = 2)
  
  # set.seed(20200814)
  # Perform analysis
  for (r in 1:reps) {
    Xy <- genData(conds[1, ], parms)
    # Xy_mis <- imposeMiss(Xy, parms)$Xy_miss
    Xy_mis <- imposeMiss_int(Xy, parms, conds[1,])
    
    O <- !is.na(Xy_mis) # matrix index of observed values
    
    # LM
    lm_GS <- lm(parms$formula, data = Xy)
    # lm_CC <- lm(parms$formula, data = Xy_mis, na.action = na.omit)
    lm_CC <- lm(parms$formula, data = Xy_mis)
    lm_sndt <- list(GS = lm_GS, CC = lm_CC)
    
    lm_par <- as.data.frame(lapply(lm_sndt, coef))
    
    R2[r, ] <- sapply(lm_sndt, function(x) summary(x)$r.squared)
    
    b_fl[r, ] <- lm_par[, "GS"]
    b_ms[r, ] <- lm_par[, "CC"]
  }
  
  # Effect of missingness on analysis
  # MCMC Estimates
  out_lm <- data.frame(full = round( colMeans(b_fl), 3),
                       miss = round( colMeans(b_ms, na.rm = TRUE), 3))
  out_lm
  
  # Bias (in terms of percentage of true value)
  BPR_lm <- round(abs(out_lm$full - out_lm$miss)/out_lm$full*100, 3)
  BPR_lm
  
# IURR Error --------------------------------------------------------------
# Run IURR many times to see if an error occurs due to varibale selection
  set.seed(1234)
  for (i in 1:1e3) {
    print(i)
    Xy <- genData(conds[1, ], parms)
    
    Xy_mis <- imposeMiss(Xy, parms)$Xy_miss
    
    imp_IURR_lasso <- impute_IURR(Xy_mis = Xy_mis, 
                                  chains = parms$chains, 
                                  iters = parms$iters,
                                  reg_type = "lasso",
                                  parms = parms)
  }

# Check data-generation ---------------------------------------------------
# Compare with same code in check of DengEtAl2016 to see that the data 
# generated is the same when using same parameters
  rm(list=ls())
  getwd()
  source("./init.R")
  set.seed(1234)
  Xy <- genData(conds[1, ], parms)
  colMeans(Xy)
  Xy_mis <- imposeMiss(Xy, parms)$Xy_miss
  mean(rowSums(is.na(Xy_mis)) != 0) # check correct miss %
  which(rowSums(is.na(Xy_mis[])) == 0)
  
  # Try imputation with DURR with same seed on same data to check 
  # equivalence of code
  imp_DURR_lasso <- impute_DURR(Xy_mis = Xy_mis,
                                chains = parms$chains,
                                iters = parms$iters,
                                reg_type="lasso",
                                parms = parms)
  imp_DURR_lasso$dats$`1`$z1
  
  imp_IURR_lasso <- impute_IURR(Xy_mis = Xy_mis,
                                chains = parms$chains,
                                iters = parms$iters,
                                reg_type = "lasso",
                                parms = parms)
  imp_IURR_lasso$dats$`1`$z1
  getwd()
  