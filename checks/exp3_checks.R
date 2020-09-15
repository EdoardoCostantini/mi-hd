### Title:    Checking parts of experiment 3 work as expected
### Project:  Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-08-31

rm(list=ls())
source("./init_general.R")
source("./exp3_init.R")

# Effects of missingness imposition on estiamtes --------------------------
# Check the effect of the missingness on the coefficients estiamtes. I see
# that the biasing of the coefficients requires the variables defining the
# prob of missingness need to be involved in the analysis model.

# > Part 1: Simple case ---------------------------------------------------

rm(list = ls())
source("./init_general.R")
source("./exp3_init.R")

# Set params
# parms$n <- 1e4 # changes nothing
# parms$lm_y_x <- 1:2 # changes nothing
# parms$int_cen <- TRUE # changes very little so keep FALSE
# parms$z_m_id <- "z1"# increases bias but does not reduce coverage
# parms$item_mean <- 0 # if not 0, CI for CC will be a bit better
# parms$rm_x <- c("y")
cond <- conds[2, ]

# Simulation Set up
set.seed(202019)
reps <- 1e4
CD_res <- CC_res <- NULL
ci.store <- vector("list", reps)

pb <- txtProgressBar(min = 0, max = reps, style = 3)

for (r in 1:reps) {
  # Gen Data
  Xy     <- simData_int(parms, cond)
  
  # Impose missing values
  Xy_mis <- imposeMiss_int(Xy, parms, cond)
  
  # Fit stuff
  fit.list <- lapply(list(Xy = Xy,
                          Xy_mis = Xy_mis), function(x) lm(y ~ z1 * z2 + z3 + z4, data = x))
  coef.list <- lm_EST(fit.list)#[-1, ]
  ci.list   <- lm_CI(fit.list)#[-c(1, 7), ]
  
  # Estiamtes
  CD_res <- cbind(CD_res, coef.list[-1,1])
  CC_res <- cbind(CC_res, coef.list[-1,2])
  
  # CI coverage
  ci.store[[r]] <- ci.list[-which(rownames(ci.list) %in% "(Intercept)"), ]
  
  # Monitoring
  setTxtProgressBar(pb, r)
}

# MCMC Estiamtes
  res <- sapply(list(CD = CD_res, 
                     CC = CC_res), 
                rowMeans)

# Confidence Interval Coverage
  ci.list <- list(CD = NULL,
                  CC = NULL)
  for (i in 1:length(ci.store)) {
    str_thrs <- nrow(ci.store[[i]])/2
    ci_low   <- ci.store[[i]][1:str_thrs, ]
    ci_hig   <- ci.store[[i]][-(1:str_thrs), ]
    for (j in 1:ncol(ci_low)) {
      ci.list[[j]] <- rbind(ci.list[[j]], ci_low[, j] < res[, 1] & res[, 1] < ci_hig[, j])
    }
  }
  
  round((res[, 2] - res[, 1])/res[, 1]*100, 0)
  round(sapply(ci.list, colMeans)*100, 0)

# > Part 2: Compare with JAV imputation -----------------------------------

  rm(list=ls())
  source("./init_general.R")
  source("./exp3_init.R")
  
  # Tweak paramters
  i     <- 4
  cond  <- conds[i, ]
  reps  <- 2
  r     <- 1
  HDrun <- FALSE
  
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
    
    pb <- txtProgressBar(min = 0, max = reps, style = 3)
    
    # Perform analysis
    for (r in 1:reps) {
      # Gen one fully-obs data
      Xy     <- simData_int(parms, cond)
      
      # Impose missing values
      Xy_mis <- imposeMiss_int(Xy, parms, cond)
      
      # Generate Interaction terms
      if(cond$p < parms$n){ 
        # LOW DIM
        interact <- computeInteract(scale(Xy_mis,
                                          center = parms$int_cen,
                                          scale  = FALSE),
                                    idVars = colnames(Xy_mis),
                                    ordVars = NULL,
                                    nomVars = NULL,
                                    moderators = colnames(Xy_mis))
        Xy_input <- cbind(Xy_mis, interact)
      } else {
        # HIGH DIM (only known interaction for optimal model)
        int    <- apply(scale(Xy_mis[, parms$zInt_id], # compute interaction term
                              center = parms$int_cen,
                              scale  = FALSE),
                        1, prod)
        Xy_input <- data.frame(Xy_mis, z1.z2 = int)
      }
      
      # Missing data
      miss_descrps <- colMeans(is.na(Xy_mis)[, parms$z_m_id])
      
      # Generate an column indexing vector based on condition
      col <- indexing_columns(Xy_input, cond)
      
      ## SINGLE DATA ##
      lm_SDin <- lapply(list(GS = Xy,
                             CC = Xy_mis), 
                        add_int_term, 
                        parms = parms
      )
      lm_sndt <- fit_lm_models(lm_SDin, mod = col$lm_mod)
      
      GS_b  <- rbind(GS_b, coef(lm_sndt$GS))
      GS_se <- rbind(GS_se, summary(lm_sndt$GS)$coefficients[, 2])
      CC_b  <- rbind(CC_b, coef(lm_sndt$CC))
      CC_se <- rbind(CC_se, summary(lm_sndt$CC)$coefficients[, 2])
      
      ## MULTIPLE IMPUTATION ##
      # High Dimensional method 
      imp_HD <- impute_DURR(Z = Xy_input[, col$CIDX],
                            O = data.frame(!is.na(Xy_input[, col$CIDX])),
                            reg_type = "lasso",
                            cond = cond,
                            perform = HDrun,
                            parms = parms)
      
      # OPTIMAL MICE Ignoring Interaction
      imp_MICE_NOV <- impute_MICE_OP(Z = Xy_mis,
                                     O = data.frame(!is.na(Xy_mis)),
                                     cond = cond,
                                     perform = TRUE,
                                     parms = parms)
      
      # OPTIMAL MICE Not Ignoring interaction
      imp_MICE_JAV <- impute_MICE_OP(Z = Xy_input[, col$CIDX_MOP],
                                     O = data.frame(!is.na(Xy_input[, col$CIDX_MOP])),
                                     cond = cond,
                                     perform = TRUE,
                                     parms = parms)
      lm_MIin <- lapply(list(imp_HD = imp_HD$dats,
                             OP_NOV = imp_MICE_NOV$dats,
                             OP_JAV = imp_MICE_JAV$dats), 
                        function(x) lapply(x, 
                                           add_int_term, 
                                           parms = parms))
      
      # Fit models
      lm_fits <- lapply(lm_MIin, fit_lm_models, mod = col$lm_mod)
      
      # Pool
      if(HDrun == TRUE){
        MI_HD_b  <- rbind(MI_HD_b, lm_pool_EST_f(lm_fits$imp_HD))
      } else {
        MI_HD_b  <- rbind(MI_HD_b, rep(NA, length(lm_pool_EST_f(lm_fits$OP_NOV))))
      }
      MI_NOV_b <- rbind(MI_NOV_b, lm_pool_EST_f(lm_fits$OP_NOV))
      MI_JAV_b <- rbind(MI_JAV_b, lm_pool_EST_f(lm_fits$OP_JAV))
      
      # Monitoring
      setTxtProgressBar(pb, r)
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
  parms$mice_ndt   <- 10
  parms$mice_iters <- 20
  
  # Use function for single condition
  set.seed(1234)
  check_missImpact(i = 2, reps = 500, 
                   conds, parms, 
                   HDrun = FALSE)

  # Parallel apply function to all conditions of interest
  out <- mclapply(X        = 1:4,
                  FUN      = check_missImpact,
                  reps     = 5e2,
                  conds    = conds,
                  parms    = parms,
                  HDrun    = FALSE,
                  mc.cores = ( 4 ) )
  
  lapply(out, function(x) round(x$BPR_lm, 1))
  lapply(out, function(x) round(x$out_lm, 3))
  out

# Convergence for single imputation ----------------------------------------
# PCA method needs needs this single imputation run before PCs are extracted.
  rm(list=ls())
  source("./init_general.R")
  source("./exp3_init.R")
  
  set.seed(1234)
  convCheck <- mclapply(X = 1:5,
                        FUN = function(i = 8){
                          cond   <- conds[8, ]
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
                        mc.cores = 5
  )

  # Save output
  store.out <- list(convCheck = convCheck,
                    parms = parms,
                    cond  = conds[8,],
                    note  = "convCheck contains a list of mids objects to check
                    that the single imputation of the target variables and interaction
                    terms obtained for exp3 DA strategy converges. 5 imputations are done
                    in parallel chans to see after how many iterations they converge. In
                    the practical use, only 1 chain (so only 1 imptuation) will be 
                    performed")
  
  saveRDS(store.out, "../checks/exp3_checks_SI_convCheck.rds") 
  
  # Read Output
  store.read <- readRDS("../checks/exp3_checks_SI_convCheck.rds")
  
  # Analyse Output
  lapply(store.read$convCheck, plot,
         y = c("z1", "z2", "z3", "z1.z2"),
         layout=c(2, 4))
  
  lapply(store.read$convCheck, plot,
         y = c("y.z1", "z2.z3", "z1.z2", "z3.z4"),
         layout=c(2, 4))
  # 200 iterations for single imputation should be fine based
  # on these results

# mcapply run -------------------------------------------------------------

  rm(list=ls())
  source("./init_general.R")
  source("./exp3_init.R")
  
  # Fix parameters to make manageble task
  parms$dt_rep     <- 2# 500 replications for averaging results
  parms$chains     <- 1 # 1   number of parallel chains for convergence check
  parms$iters      <- 2 # 75  total iterations
  parms$burnin_imp <- 0 # 50  how many imputation iterations discarded
  parms$ndt        <- 2 # 10  number of imputed datasets to pool esitmaes from
  parms$thin       <- (parms$iters - parms$burnin_imp)/parms$ndt
  parms$pos_dt  <- (parms$burnin_imp+1):parms$iters # candidates
  parms$keep_dt <- parms$pos_dt[seq(1, 
                                    length(parms$pos_dt), 
                                    parms$thin)] # keep 1 dataset every thin
  parms$chains_bl     <- 1 # 1 
  parms$iters_bl      <- 2 # 2e3  total iterations
  parms$burnin_imp_bl <- 0 # 1950 discarded iterations
  parms$thin_bl       <- (parms$iters_bl - parms$burnin_imp_bl)/parms$ndt
  parms$pos_dt_bl     <- (parms$burnin_imp_bl+1):parms$iters_bl # candidate
  parms$keep_dt_bl    <- parms$pos_dt_bl[seq(1, 
                                             length(parms$pos_dt_bl), 
                                             parms$thin_bl)]
  parms$mice_iters <- 2 #  20
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
                   imp_values   = TRUE)
  
  # Single run of doRep
  check_doRep <- doRep(1, conds = conds, parms = parms, debug = FALSE)
  
  lapply(check_doRep, 
         function(x) round(x$lm_EST, 1)) # everything is there
  lapply(check_doRep[5:8], 
         function(x) round(x$lm_EST, 1)) # DA includes the two version of DA
  lapply(check_doRep[-(5:8)], 
         function(x) round(x$lm_EST, 1)) # not DA have the regular methods
  
  # Run mcapply
  out <- mclapply(X        = 1 : parms$dt_rep,
                  FUN      = doRep,
                  conds    = conds,
                  parms    = parms,
                  debug    = FALSE,
                  mc.cores = ( parms$dt_rep ) )
  out$parms <- parms
  out$conds <- conds
  
  # Check results presence
  # Time
  out_time <- sapply((1:nrow(out$conds))[-(5:8)],
                     res_sem_time, out = out)
  colnames(out_time) <- paste0("cond_", 1:nrow(out$conds))[-(5:8)]
  
  out_time_DA <- sapply((1:nrow(out$conds))[(5:8)], 
                     res_sem_time, out = out)
  colnames(out_time_DA) <- paste0("cond_", 1:nrow(out$conds))[(5:8)]
  
  out_time
  out_time_DA
  
  ## SEM estiamtes raw data (saturated model) ##
  # Extract results per conditions
  sem_res <- lapply(1:nrow(out$conds),
                     function(x) res_sum(out, 
                                         model = "sem", 
                                         condition = x))
  
  lm_res <- lapply(1:nrow(out$conds),
                    function(x) res_sum(out, 
                                        model = "lm", 
                                        condition = x))
  
  # Show results all conditions
  lapply(1:nrow(out$conds),
         function(x) sem_res[[x]]$bias_per)
  
  lapply(1:nrow(out$conds),
         function(x) lm_res[[x]]$bias_per)
  
  