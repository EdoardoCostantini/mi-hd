### Title:    Replication Zhao Long 2016 - Check certain key elements
### Author:   Edoardo Costantini
### Created:  2020-02-20
### Modified: 2020-02-20
### Notes:    Check whether some setups of the simulation are achieved

source("./init.R")
source("./functions.R")

# Missingness -------------------------------------------------------------

  cond <- conds[1,]
  obj <- NULL
  for (i in 1:1e3) {
    Xy <- genData(cond, parms)
    obj[i] <- mean(imposeMiss(Xy, parms)$nR)
  }
  mean(obj)
  # Should be .3 missingness on this variable on average
  

# GS and CC check ---------------------------------------------------------
# results of table 1 that do not require imputation
  reps <- 500
  bias_GS <- NULL
  CI_GS <- matrix(rep(NA, 4*reps), ncol = 4)
  bias_CC <- NULL
  CI_CC <- matrix(rep(NA, 4*reps), ncol = 4)
  
  for (i in 1:reps) {
    Xy <- genData(conds[1, ], parms)
    mod_GS <- lm(parms$formula, data = Xy)
    bias_GS[i] <- mod_GS$coefficients["z1"] - 1
    CI_GS[i, ] <- check_cover(confint(mod_GS))
    
    Xy_mis <- imposeMiss(Xy, parms)$Xy_miss
    mod_CC <- lm(parms$formula, data = Xy_mis)
    bias_CC[i] <- mod_CC$coefficients["z1"] - 1
    CI_CC[i, ] <- check_cover(confint(mod_CC))
  }
  round(c(mean(bias_GS), mean(bias_CC)), 3)
  round(c(sd(bias_GS), sd(bias_CC)), 3)
  colMeans(CI_GS)
  colMeans(CI_CC)
  
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
  
  # Try imputation with DURR with same seed on same data to check 
  # equivalence of code
  imp_DURR_lasso <- impute_DURR(Xy_mis = Xy_mis,
                                chains = parms$chains,
                                iters = parms$iters,
                                reg_type="lasso")
  imp_DURR_lasso$`1`$z1
  
  imp_IURR_lasso <- impute_IURR(Xy_mis = Xy_mis, 
                                cond = cond, 
                                chains = parms$chains, 
                                reg_type = "lasso")
  imp_IURR_lasso$`1`$z1
  getwd()