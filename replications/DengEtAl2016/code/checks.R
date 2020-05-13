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
    
    # # TRUE imputation model used
    # imp_MICE_TR <- impute_MICE_TR(Xy_mis = Xy_mis, 
    #                               AS_size = conds[1, 3],
    #                               chains = parms$chains, 
    #                               iters = parms$iters, 
    #                               parms = parms)
    # fits_md <- fit_models(imp_MICE_TR$dats, parms$formula)
    # # Same as
    # # fits_md <- mice::lm.mids(parms$formula, data = imp_MICE_TR$mids)
    # # poolFit1 <- mice::pool(fits1)
    # # summary(poolFit1)[, 1]
    # bias_TR[i, ] <- get_pool_EST(fits_md) - 1
    # CI_TR[i, ] <- check_cover(get_pool_CI(fits_md))
  }
  # Average bias
  round(
    t(data.frame(GS = colMeans(bias_GS),
                 CC = colMeans(bias_CC),
                 TR = colMeans(bias_TR)
                 )), 3)
  
  # Average Standard Error
  round(
    t(data.frame(GS = colMeans(se_GS),
                 CC = colMeans(se_CC)
                 # TR = colMeans(se_TR)
    )), 3)
  
  # Monte Carlo SD
  round(
    t(data.frame(GS = apply(coef_GS, 2, sd),
                 CC = apply(coef_CC, 2, sd)
    )), 3)
  
  # Average Coverage
  t(sapply(list(CI_GS,CI_CC, CI_TR), colMeans))

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
  set.seed(1234)
  Xy <- genData(conds[1, ], parms)
  colMeans(Xy)
  Xy_miss <- imposeMiss(Xy, parms)$Xy_miss
  mean(rowSums(is.na(Xy_miss)) != 0) # check correct miss %
  which(rowSums(is.na(Xy_mis[])) == 0)
