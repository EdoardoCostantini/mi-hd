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
  reps <- 10
  r <- 1
  # cond$r2 <- .9
  # cond$pm <- .7
  # parms$missType <- "tails"
  
  # Make it a function
  check_missImpact <- function(i, reps, conds, parms) {
    
    # For replicability
    .lec.SetPackageSeed(rep(parms$seed, 6))
    if(!1 %in% .lec.GetStreams()) # if the streams do not exist yet
      .lec.CreateStream(c(1 : parms$nStreams))
    .lec.CurrentStream(1) # this is equivalent to setting the seed
    
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
        col_int <- paste0("z", parms$yMod_int)
        
        # Fully obsrved
        int_term <- apply(scale(Xy[, col_int], 
                                center = TRUE,
                                scale = FALSE), 
                          1, 
                          prod)
        
        Xy <- cbind(Xy, int_term)
        colnames(Xy)[colnames(Xy) == "int_term"] <- paste0(col_int, 
                                                           collapse = "")
        
        # W/ missings
        int_term <- apply(scale(Xy_mis[, col_int], 
                                center = TRUE,
                                scale = FALSE), 
                          1, prod)
        Xy_mis <- cbind(Xy_mis, int_term)
        colnames(Xy_mis)[colnames(Xy_mis) == "int_term"] <- paste0(col_int, 
                                                                   collapse = "")
      }
      
      # Miss description
      O <- is.na(data.frame(Xy_mis))
      
      if(parms$zm_n > 1){
        pm[r,]<- colMeans(O[, (1:parms$zm_n)+1])
      } else {
        pm[r,] <- mean(O[, (1:parms$zm_n)+1])
      }

      # Fit models and store output
      # Multiple datasets
      if(cond$int_sub == FALSE){
        mod <- parms$frm
      } else {
        mod <- paste0("y ~ -1 + ",
                      paste0("z", parms$yMod_cov, collapse = " + "),
                      " + ",
                      paste0("z", parms$yMod_int, collapse = ""))
      }
      lm_sndt <- fit_lm_models(list(GS      = Xy,
                                    CC      = Xy_mis),
                               mod = mod)
      
      GS_b  <- rbind(GS_b, coef(lm_sndt$GS))
      GS_se <- rbind(GS_se, summary(lm_sndt$GS)$coefficients[, 2])
      CC_b  <- rbind(CC_b, coef(lm_sndt$CC))
      CC_se <- rbind(CC_se, summary(lm_sndt$CC)$coefficients[, 2])
      
      # Imputation
      imp_MICE_OP <- impute_MICE_OP(Z = Xy_mis,
                                    O = O,
                                    cond = cond,
                                    perform = parms$meth_sel$MI_OP,
                                    parms = parms)
      lm_fits <- fit_lm_models(imp_MICE_OP$dats , mod = mod)
      lm_pool_est <- lm_pool_EST_f(lm_fits)
      MI_b <- rbind(MI_b, lm_pool_est)
      
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
    return(list(BPR_lm=BPR_lm,
                out_lm=out_lm))
  }
  
  # Use function for single condition
  check_missImpact(i = 2, reps = 1e2, conds, parms)
  
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


# mcapply run -------------------------------------------------------------

  rm(list=ls())
  source("./init_general.R")
  source("./exp3_init.R")
  
  # Fix parameters to manageble task
  parms$dt_rep     <- 5# 500 replications for averaging results
  parms$chains     <- 1 # 1   number of parallel chains for convergence check
  parms$iters      <- 5 # 75  total iterations
  parms$burnin_imp <- 0 # 50  how many imputation iterations discarded
  parms$ndt        <- 5 # 10  number of imputed datasets to pool esitmaes from
  parms$thin       <- (parms$iters - parms$burnin_imp)/parms$ndt
  parms$pos_dt  <- (parms$burnin_imp+1):parms$iters # candidates
  parms$keep_dt <- parms$pos_dt[seq(1, 
                                    length(parms$pos_dt), 
                                    parms$thin)] # keep 1 dataset every thin
  parms$chains_bl     <- 1 # 1 
  parms$iters_bl      <- 5 # 2e3  total iterations
  parms$burnin_imp_bl <- 0 # 1950 discarded iterations
  parms$thin_bl       <- (parms$iters_bl - parms$burnin_imp_bl)/parms$ndt
  parms$pos_dt_bl     <- (parms$burnin_imp_bl+1):parms$iters_bl # candidate
  parms$keep_dt_bl    <- parms$pos_dt_bl[seq(1, 
                                             length(parms$pos_dt_bl), 
                                             parms$thin_bl)]
  parms$mice_iters <- 5 #  20
  parms$mice_ndt   <- parms$ndt
  
  parms$store <- c(cond         = TRUE,
                   dat_full     = FALSE,
                   dat_miss     = FALSE,
                   sem_EST      = TRUE,
                   sem_CI       = TRUE,
                   lm_EST       = TRUE,
                   lm_CI        = TRUE,
                   miss_descrps = TRUE,
                   run_time_min = TRUE,
                   imp_values   = FALSE)
  
  # Run mcapply
  
  out <- mclapply(X        = 1 : parms$dt_rep,
                  FUN      = doRep,
                  conds    = conds,
                  parms    = parms,
                  debug    = FALSE,
                  mc.cores = ( parms$dt_rep ) )
  out$parms <- parms
  
  # Apply some rsults functions
  
  condition <- 1
  select_cond <- names(out[[1]])[condition]
  
  # Time
  out_time <- sapply(1:length(names(out[[1]])), res_sem_time, out = out)
  colnames(out_time) <- names(out[[1]])
  t(out_time)
  
  # Estaimtes Analysis ------------------------------------------------------
  out[[1]]$`cond_0.3_1e-05_0.65_FALSE_FALSE_FALSE_25`$sem_EST
  out[[1]]$`cond_0.3_1e-05_0.65_FALSE_FALSE_FALSE_25`$lm_EST
  ## SEM estiamtes raw data (saturated model) ##
  # Extract results per conditions
  sem_res <- lapply(1:length(out[[1]]),
                     function(x) res_sum(out, 
                                         model = "sem", 
                                         condition = x))
  
  lm_res <- lapply(1:length(out[[1]]),
                    function(x) res_sum(out, 
                                        model = "lm", 
                                        condition = x))
  
  # Show results all conditions for a given data rep
  lapply(1:length(out[[1]]),
         function(x) sem_res[[x]]$bias_per)
  
  lapply(1:length(out[[1]]),
         function(x) lm_res[[x]]$bias_per)
  
  