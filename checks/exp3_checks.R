### Title:    Checking parts of experiment 3 work as expected
### Project:  Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-05-19

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

# > explaied variance on averegae -------------------------------------------

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
  
# > effects of missingness imposition on estiamtes --------------------------
  
  rm(list=ls())
  source("./init_general.R")
  source("./exp3_init.R")
  
  # Tweak paramters
  i     <- 2
  cond  <- conds[i, ]
  reps  <- 2
  r     <- 1
  HDrun <- FALSE
  
  # cond$r2 <- .9
  # cond$pm <- .7
  set.seed(1234)
  
  # Make it a function
  check_missImpact <- function(i, reps, conds, parms, HDrun = TRUE) {
    
    # For replicability
    .lec.SetPackageSeed(rep(parms$seed, 6))
    if(!1 %in% .lec.GetStreams()) # if the streams do not exist yet
      .lec.CreateStream(c(1 : parms$nStreams))
    .lec.CurrentStream(1) # this is equivalent to setting the seed
    
    # Select condition
    cond <- conds[i, ]

    # Create/empty storing objects
    GS_b <- CC_b <- GS_se <- CC_se <- MI_b <- NULL
    MI_NOV_b <- MI_JAV_b <- MI_HD_b <- NULL
    pm <- matrix(NA, nrow = reps, ncol = parms$zm_n)
    
    # Perform analysis
    for (r in 1:reps) {
      Xy     <- simData_int(parms, cond)
      Xy_mis <- imposeMiss_int(Xy, parms, cond)
      
      # Cast data to method specific format
      # For Optimal Impuation
      col_int  <- paste0("z", parms$yMod_int)    # variables interacting
      z1.z2    <- apply(scale(Xy_mis[, col_int], # compute interaction term
                              center = FALSE,
                              scale = FALSE),
                        1, prod)
      Xy_MIOP <- cbind(Xy_mis, z1.z2)
      
      # For all methods
      # Generate Interaction terms (for low dimensional condition)
      if(cond$p < parms$n){
        interact <- computeInteract(Xy_mis,
                                    idVars = colnames(Xy_mis),
                                    ordVars = NULL,
                                    nomVars = NULL,
                                    moderators = colnames(Xy_mis))
        Xy_input <- cbind(Xy_mis, interact)
      } else {
        Xy_input <- Xy_mis
      }
      
      # Populate it with values if there are missings (single imputation)
      if(cond$int_da == TRUE){
        predMatrix <- quickpred(Xy_input, mincor = .3)
        Xy_mis_DAmids  <- mice(Xy_input,
                               m               = 1, 
                               maxit           = parms$SI_iter,
                               predictorMatrix = predMatrix,
                               # printFlag       = FALSE,
                               ridge           = cond$ridge,
                               method          = "pmm")
        Xy_input <- complete(Xy_mis_DAmids)
        Xy_input[, which(names(Xy) %in% parms$z_m_id)] <- 
          Xy_mis[, which(names(Xy) %in% parms$z_m_id)]
      }
      
      # Missing data 
      # O <- !is.na(Xy_mis) # matrix index of observed values
      miss_descrps <- colMeans(is.na(Xy_mis)[, parms$z_m_id])
      
      # Generate an column indexing vector based on condition
      if(cond$int_sub == FALSE & cond$int_da == FALSE){
        CIDX_all <- colnames(Xy_input)[!grepl("\\.", 
                                              colnames(Xy_input))]
        CIDX_MOP <- colnames(Xy_MIOP)[!grepl("\\.", 
                                             colnames(Xy_MIOP))]
        mod_MICE_OP <- mod <- parms$frm
      }
      if(cond$int_sub == TRUE & cond$int_da == FALSE){
        CIDX_all <- colnames(Xy_input)[!grepl("\\.", 
                                              colnames(Xy_input))]
        CIDX_MOP <- c(colnames(Xy_MIOP)[!grepl("\\.", 
                                               colnames(Xy_MIOP))],
                      "z1.z2")
        mod         <- parms$frm_int
        mod_MICE_OP <- paste0("y ~ -1 + ",
                              paste0("z", parms$lm_y_x, collapse = " + "),
                              " + ",
                              paste0("z", parms$lm_y_i, collapse = "."))
      }
      if(cond$int_sub == FALSE & cond$int_da == TRUE){
        CIDX_all <- colnames(Xy_input)
        CIDX_MOP <- colnames(Xy_MIOP)[!grepl("\\.", 
                                             colnames(Xy_MIOP))]
        mod_MICE_OP <- mod <- parms$frm
      }
      if(cond$int_sub == TRUE & cond$int_da == TRUE){
        CIDX_all <- colnames(Xy_input)
        CIDX_MOP <- c(colnames(Xy_MIOP)[!grepl("\\.", 
                                               colnames(Xy_MIOP))],
                      "z1.z2")
        mod         <- parms$frm_int
        mod_MICE_OP <- paste0("y ~ -1 + ",
                              paste0("z", parms$lm_y_x, collapse = " + "),
                              " + ",
                              paste0("z", parms$lm_y_i, collapse = "."))
      }
      
      ## SINGLE DATA ##
      lm_sndt <- fit_lm_models(list(GS      = Xy,
                                    CC      = Xy_mis),
                               mod = mod)
      
      GS_b  <- rbind(GS_b, coef(lm_sndt$GS))
      GS_se <- rbind(GS_se, summary(lm_sndt$GS)$coefficients[, 2])
      CC_b  <- rbind(CC_b, coef(lm_sndt$CC))
      CC_se <- rbind(CC_se, summary(lm_sndt$CC)$coefficients[, 2])
      
      ## MULTIPLE IMPUTATION ##
      # High Dimensional method 
      imp_HD <- impute_DURR(Z = Xy_input[, CIDX_all],
                            O = data.frame(!is.na(Xy_input[, CIDX_all])),
                            reg_type = "lasso",
                            cond = cond,
                            perform = HDrun,
                            parms = parms)
      
      # OPTIMAL MICE Ignoring Interaction
      imp_MICE_OP_NOV <- impute_MICE_OP(Z       = Xy_mis,
                                        O       = data.frame(!is.na(Xy_mis)),
                                        cond    = cond,
                                        perform = TRUE,
                                        parms   = parms)
      
      # OPTIMAL MICE Not Ignoring interaction
      imp_MICE_OP_JAV <- impute_MICE_OP(Z = Xy_MIOP[, CIDX_MOP],
                                        O = data.frame(!is.na(Xy_MIOP[, CIDX_MOP])),
                                        cond = cond,
                                        perform = TRUE,
                                        parms = parms)
      
      # Fit Models
      MICE_HD <- fit_lm_models(imp_HD$dats,
                               mod = mod)
      MICE_OP_NOV <- fit_lm_models(imp_MICE_OP_NOV$dats,
                                   mod = mod)
      MICE_OP_JAV <- fit_lm_models(imp_MICE_OP_JAV$dats,
                                   mod = mod_MICE_OP)
      
      # Pool
      if(HDrun == TRUE){
        MI_HD_b  <- rbind(MI_HD_b, lm_pool_EST_f(MICE_HD))
      } else {
        MI_HD_b  <- rbind(MI_HD_b, rep(NA, length(lm_pool_EST_f(MICE_OP_NOV))))
      }
      MI_NOV_b <- rbind(MI_NOV_b, lm_pool_EST_f(MICE_OP_NOV))
      MI_JAV_b <- rbind(MI_JAV_b, lm_pool_EST_f(MICE_OP_JAV))
    }
    
    # Effect of missingness on analysis
    # MCMC Estimates
    out_lm <- t(data.frame(GS = round( colMeans(GS_b), 3),
                           CC = round( colMeans(CC_b, na.rm = TRUE), 3),
                           NOV = round( colMeans(MI_NOV_b, na.rm = TRUE), 3),
                           JAV = round( colMeans(MI_JAV_b, na.rm = TRUE), 3),
                           HD = round( colMeans(MI_HD_b, na.rm = TRUE), 3) ))
    out_lm
    
    # Bias (in terms of percentage of true value)
    BPR_lm <- t(data.frame(CC  = (out_lm["GS",] - out_lm["CC",])/out_lm["GS",] * 100,
                           NOV = (out_lm["GS",] - out_lm["NOV",])/out_lm["GS",] * 100,
                           JAV = (out_lm["GS",] - out_lm["JAV",])/out_lm["GS",] * 100,
                           HD  = (out_lm["GS",] - out_lm["HD",])/out_lm["GS",] * 100))
    BPR_lm
    
    # Store it
    return(list(out_lm = out_lm,
                BPR_lm = BPR_lm))
  }
  
  # Modify parameters for check
  parms$iters      <- 75 # 75  total iterations
  parms$burnin_imp <- 50 # 50  how many imputation iterations discarded
  parms$ndt        <- 10 # 10  number of imputed datasets to pool esitmaes from
  parms$thin       <- (parms$iters - parms$burnin_imp)/parms$ndt
  parms$pos_dt     <- (parms$burnin_imp+1):parms$iters # candidates
  parms$keep_dt    <- parms$pos_dt[seq(1, 
                                       length(parms$pos_dt), 
                                       parms$thin)] # keep 1 dataset every thin
  parms$n          <- 1e3
  parms$mice_ndt   <- 10
  parms$mice_iters <- 20
  parms$SI_iter    <- 3e2
  parms$missType <- c("high", "low", "tails")[1]
  
  # Use function for single condition
  check_missImpact(i = 1, reps = 50, conds, parms, HDrun = FALSE)
  
  # Parallel apply function to all conditions of interest
  # parms$missType <- "tails"
  out <- mclapply(X        = 5:8,
                  FUN      = check_missImpact,
                  reps     = 50,
                  conds    = conds,
                  parms    = parms,
                  HDrun    = TRUE,
                  mc.cores = ( 4 ) )
  lapply(out, function(x) round(x$BPR_lm, 1))
  lapply(out, function(x) round(x$out_lm, 3))
  out


# Convergence for single imputation ----------------------------------------
  
  rm(list=ls())
  source("./init_general.R")
  source("./exp3_init.R")
  
  set.seed(1234)
  convCheck <- mclapply(X = 1:3,
                        FUN = function(i = 8){
                          cond <- conds[8, ]
                          Xy     <- simData_int(parms, cond)
                          Xy_mis <- imposeMiss_int(Xy, parms, cond)
                          
                          # For all methods
                          # Generate Interaction terms (for low dimensional condition)
                          interact <- computeInteract(Xy_mis,
                                                      idVars = colnames(Xy_mis),
                                                      ordVars = NULL,
                                                      nomVars = NULL,
                                                      moderators = colnames(Xy_mis))
                          Xy_input <- cbind(Xy_mis, interact)
                          
                          # Populate it with values if there are missings (single imputation)
                          predMatrix <- quickpred(Xy_input, mincor = .3)
                          Xy_mis_DAmids  <- mice(Xy_input,
                                                 m               = 5,
                                                 maxit           = 500,
                                                 predictorMatrix = predMatrix,
                                                 printFlag       = FALSE,
                                                 ridge           = cond$ridge,
                                                 method          = "norm")
                          return(Xy_mis_DAmids)
                        },
                        mc.cores = 3
  )
  lapply(convCheck, plot,
         y = c("z1", "z2", "z1.z2", "z3.z4"),
         layout=c(2, 4))
  
  lapply(convCheck, plot,
         y = c("y.z1", "z2.z3", "z1.z2", "z3.z4"),
         layout=c(2, 4))
  # Settlede on a 300 iterations for single imputation

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
  
  # Analysis
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
  
  