### Title:    Checking parts of experiment 3 work as expected
### Project:  Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-05-19

# Data Generation Latent Structure ----------------------------------------

  rm(list=ls())
  source("./init_general.R")
  source("./exp3_init.R")

# > works for all conditions ----------------------------------------------

  all_conds_data <- lapply(1:nrow(conds), function(x){
    X <- simData_lv(parms, conds[x,])
    Xy_mis <- imposeMiss_lv(X, parms, conds[x, ])
    return(Xy_mis)
  })
  length(all_conds_data)

# Explaied variance on averegae -------------------------------------------

  cond <- conds[1, ]
  res <- NULL
  
  
  for (i in 1:1e3) {
    Xy <- simData_int(parms, cond)
    
    if(cond$int_sub == FALSE){
      lm_fit <- summary(lm(parms$frm, data = Xy))
    } else {
      lm_fit <- summary(lm(parms$frm_int, data = Xy))
    }
    res[i] <-  lm_fit$r.squared
  }
  
  c(true = cond$r2, 
    avg = mean(res))
  
# Effects of missingness imposition on estiamtes --------------------------
  
  rm(list=ls())
  source("./init_general.R")
  source("./exp3_init.R")
  
  # Tweak paramters
  cond <- conds[2, ]
  # cond$r2 <- .9
  # cond$pm <- .7
  # parms$missType <- "tails"
  
  # Make it a function
  check_missImpact <- function(i, reps, conds, parms) {
    
    # For replicability
    .lec.SetPackageSeed(rep(parms$seed, 6))
    if(!rp %in% .lec.GetStreams()) # if the streams do not exist yet
      .lec.CreateStream(c(1 : parms$nStreams))
    .lec.CurrentStream(rp) # this is equivalent to setting the seed
    
    # Progres bar
    pb = txtProgressBar(min = 0, max = reps, initial = 1,
                        title = "Check missing imposition effects",
                        style = 3)
    
    # Select condition
    cond <- conds[i, ]
    
    # Create/empty storing objects
    GS_b <- CC_b <- GS_se <- CC_se <- MI_b <- NULL
    pm <- matrix(NA, nrow = reps, ncol = parms$zm_n)
    
    # Perform analysis
    for (r in 1:reps) {
      
      # Gen data
      Xy <- simData_int(parms, cond)
      Xy_mis <- imposeMiss_int(Xy, parms, cond)
      Z <- as.matrix(Xy[, -1])
      y <- Xy$y
      
      if(cond$int_sub == TRUE){
        col_int  <- paste0("z", parms$yMod_int)
        int_term <- apply(scale(Xy_mis[, col_int], TRUE, FALSE), 1, prod)
        Xy_mis   <- cbind(Xy_mis, int_term)
        colnames(Xy_mis)[colnames(Xy_mis) == "int_term"] <- paste0(col_int, 
                                                                   collapse = "")
      }
      
      # Miss description
      O <- sapply(data.frame(Xy_mis), is.na)
      if(parms$zm_n > 1){
        pm[r,]<- colMeans(O[, (1:parms$zm_n)+1])
      } else {
        pm[r,] <- mean(O[, (1:parms$zm_n)+1])
      }
      ind_cc <- rowSums(O) == 0
      
      col_indx <- paste0("z", parms$yMod_cov)
      
      # Fit models and store output
      if(cond$int_sub == FALSE){
        lm_GS <- lm(parms$frm, data = Xy)
        lm_CC <- lm(parms$frm, data = Xy_mis)
      } else {
        lm_GS <- lm(parms$frm_int, data = Xy)
        lm_CC <- lm(parms$frm_int, data = Xy_mis)
      }
      GS_b  <- rbind(GS_b, coef(lm_GS))
      GS_se <- rbind(GS_se, summary(lm_GS)$coefficients[, 2])
      CC_b  <- rbind(CC_b, coef(lm_CC))
      CC_se <- rbind(CC_se, summary(lm_CC)$coefficients[, 2])
      
      # Impute and analyze with mice optimal
      predMat <- matrix(rep(0, ncol(Xy_mis)^2), ncol = ncol(Xy_mis), 
                        dimnames = list(colnames(Xy_mis), colnames(Xy_mis)))
      if(cond$int_sub == FALSE){
        # predMat[parms$z_m_id, c("y", "z6", "z7")] <- 1
        predMat[parms$z_m_id, 
                c(y=1, parms$blck1+1, parms$blck2+1)] <- 1
      } else {
        predMat[c(parms$z_m_id, paste0("z", parms$yMod_int, collapse = "")), 
                c(y=1, parms$blck1+1, parms$blck2+1, int=ncol(predMat))] <- 1
      }
      diag(predMat) <- 0
      
      mice_out <- mice::mice(Xy_mis,
                             m = 5,
                             maxit = 10,
                             predictorMatrix = predMat,
                             printFlag = FALSE,
                             method = "pmm")
      
      # Store results
      mice_dat <- mice::complete(mice_out, "all")
      if(cond$int_sub == FALSE){
        mice_fit <- lapply(mice_dat, lm, formula = parms$frm)
      } else {
        mice_fit <- lapply(mice_dat, 
                           lm, 
                           formula = paste0("y ~ -1 + ",
                                            paste0("z", parms$yMod_cov, collapse = " + "),
                                            " + ",
                                            paste0("z", parms$yMod_int, collapse = ""))
        )
      }
      mice_sum <- lapply(mice_fit, summary)
      mice_est <- t(sapply(mice_sum, function(x) coef(x)[, "Estimate"]))
      mice_pool <- colMeans(mice_est)
      MI_b  <- rbind(MI_b, mice_pool)
      
      setTxtProgressBar(pb, r)
    }
    # Effect of missingness on analysis
    # MCMC Estimates
    out_lm <- t(data.frame(GS = round( colMeans(GS_b), 3),
                           CC = round( colMeans(CC_b, na.rm = TRUE), 3),
                           MI = round( colMeans(MI_b), 3)))
    out_lm
    
    # Bias (in terms of percentage of true value)
    BPR_lm <- t(data.frame(CC = abs(out_lm["GS",] - out_lm["CC",])/out_lm["GS",] * 100,
                           MI = abs(out_lm["GS",] - out_lm["MI",])/out_lm["GS",] * 100))
    BPR_lm
    
    # Store it
    return(BPR_lm)
  }
  
  # Use function for single condition
  check_missImpact(i = 4, reps = 5, conds, parms)
  
  # Parallel apply function to all conditions of interest
  out <- mclapply(X        = 1 : 4,
                  # Only the first 4 conditions are important to check.
                  # The others do no change the data generation but the 
                  # treatment of the data
                  FUN      = check_missImpact,
                  reps     = 5e2,
                  conds    = conds,
                  parms    = parms,
                  mc.cores = ( 4 ) )
  