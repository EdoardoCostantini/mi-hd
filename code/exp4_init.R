# Title:    Initialization scirpt parameter
# Project:  Imputing High Dimensional Data
# Author:   Edoardo Costantini
# Created:  2020-07-01
# Modified: 2022-02-25

# Fixed Parameters --------------------------------------------------------

  parms <- list()
  parms$exp <- 4 # which experiment is been run

  # Itereations, repetitions, etc
  parms$dt_rep     <- 5 # 500 replications for averaging results
  parms$chains     <- 1  # 1   number of parallel chains for convergence check
  parms$iters      <- 10 # 70  total iterations
  parms$burnin_imp <- 5  # 60  how many imputation iterations discarded
  parms$ndt        <- 5  # 10  number of imputed datasets to pool esitmaes from
  parms$thin       <- (parms$iters - parms$burnin_imp)/parms$ndt
  parms$pos_dt     <- (parms$burnin_imp+1):parms$iters # candidates
  parms$keep_dt    <- parms$pos_dt[seq(1, 
                                       length(parms$pos_dt), 
                                       parms$thin)] # keep 1 dataset every thin
  
  # For blasso
  parms$chains_bl     <- 1  # 1 
  parms$iters_bl      <- 10 # 70 total iterations
  parms$burnin_imp_bl <- 5  # 60 discarded iterations
  parms$thin_bl       <- (parms$iters_bl - parms$burnin_imp_bl)/parms$ndt
  parms$pos_dt_bl     <- (parms$burnin_imp_bl+1):parms$iters_bl # candidate
  parms$keep_dt_bl    <- parms$pos_dt_bl[seq(1, 
                                             length(parms$pos_dt_bl), 
                                             parms$thin_bl)]
  
  # For mice-like algorithms
  parms$mice_iters <- 10 # 60
  parms$mice_ndt   <- parms$ndt # mice keeps data in a different way
  
# Data gen ----------------------------------------------------------------
  
  # "True" values
  parms$item_mean <- 0
  parms$item_var  <- 1
  
  # Variable blocks
  parms$blcky <- 1:1 # dependent variables
  
  parms$blck1 <- 1:10 # higly correlated
  parms$blck2 <- (tail(parms$blck1,1)+1):(tail(parms$blck1,1)+10) # correlated
  # y is excluded from this naming convention. Column 1 to 10 of the
  # X matrix, y is appended afterwards
  
  parms$blck1_r <- .6 # correlation for highly correlated variables
  parms$blck2_r <- .3 # correlation for correlated variables
  
  # Fully observed variables
  parms$Z_o_mu  <- 0 # mean for gen of fully observed variables
  
  # Variables that will have missingness
  parms$z_m_id  <- paste0("v", c(185:187, "174_LR", 118, 156))
  parms$zm_n    <- length(parms$z_m_id)
  
  # y gen / imporntant predictors
  parms$yMod_cov <- parms$blck1[1:4]
  parms$yMod_int <- parms$yMod_cov[1:2]
  parms$zInt_id  <- paste0("z", parms$yMod_int)
  parms$int_cen  <- FALSE
  parms$b_main   <- rep(1, length(parms$yMod_cov))
  parms$b_int    <- 1
  
  parms$lm_y_x  <- parms$blck1[1:4] # predictors for the substantive linear model
  parms$lm_y_i  <- parms$lm_y_x[1:2]
  
  parms$b       <- 1 # for every regression coefficient
  
# Models ------------------------------------------------------------------
  parms$frm     <- paste0("y ~ ",
                          paste0("z", parms$yMod_cov, collapse = " + "))
  parms$frm_int <- paste0("y ~ ",
                          paste0("z", parms$yMod_cov, collapse = " + "),
                          " + ",
                          paste0("z", parms$yMod_int, collapse = "."))

# Imputation methods ------------------------------------------------------
  parms$alphaCI <- .95 # confidence level for parameters CI
  parms$meth_sel <- data.frame(DURR_la    = FALSE,
                               IURR_la    = FALSE,
                               blasso     = FALSE,
                               bridge     = FALSE,
                               MI_PCA     = FALSE,
                               MI_CART    = FALSE,
                               MI_RF      = FALSE,
                               MI_qp      = TRUE,
                               MI_am      = TRUE,
                               MI_OP      = TRUE,
                               missFor    = TRUE,
                               mean       = TRUE,
                               CC         = TRUE,
                               GS         = TRUE)
  
  parms$methods <- names(parms$meth_sel)[which(parms$meth_sel==TRUE)]
    # (missFor, mean, CC, GS always last, alwyas present)
  
  # Location
  parms$missType <- c("high", "low", "tails")[2]
  
  # Response Model (rm)
  parms$rm_x <- c("age",          # for item- not unit-nonresponse!
                  "v243_ISCED_1", # education
                  "v35")          # trust new person
  
  # weighting the importance of predictors: all the same
  parms$auxWts <- rep(1, length(parms$rm_x))
  
  # IURR
  parms$k_IURR  <- 0 # k value to bias coef sampling covariance matrix
                     # procedure to solve possible issues of singularity

  # PCA
  parms$SI_iter      <- 10L  # 200L # iterations for single imputation in PCA run
  parms$PCA_pcthresh <- .5 # proportion of vairance for selection of PCs
  
  # Random Forest
  parms$rfntree <- 10
  parms$missFor_maxiter <- 20 # maxiter = 20
  parms$missFor_ntree <- 100  # ntree = 100
  
  # MICE true
  parms$S_all <- c(parms$rm_x,   # variables influencing the missingness
                   parms$z_m_id, # imputation of v118 will not use v118 even if its here
                   paste0("v", c(31, 126, 120:121, 127, 131, 225, 6,
                                 145, 110, 97:101, 234, 54)),
                   "v51v52_comb", "v246_egp", "v276_r")

  # Analysis model variables
  parms$am_vars <- c("v156", # attitudes toward euthanasia
                     "v31", # General Trust
                     "v126", # Confidence in Health care sys
                     "v118", # Confidence in press
                     c("v121", "v120", "v127", "v131"), # items for confidence in state scale
                     "v174_LR", # Left / Right voting
                     "country", # Country
                     "v225", # binary gender
                     "v246_egp", # socio economic status
                     c("v185", "v186", "v187"), # items making up Native attitudes
                     "v145", # attitude toward strong leaders
                     "v110", # attitude toward "order" in law and order
                     "v97", # Political Interest
                     paste0("v", 98:101), # items for politicla action scale
                     "age", # age
                     "v243_ISCED_1", # education
                     "v234", # marital status
                     "v276_r", # urb
                     "v6", # Religiousness
                     "v51v52_comb"  # Denomination
  )

# Simulation desing -------------------------------------------------------
  # Replicability
  parms$seed     <- 20220131
  parms$nStreams <- 1000

  # Output and Progres report related
  parms$outDir <- "../output/"
  parms$start_time <- format(Sys.time(), "%Y%m%d_%H%M")
  parms$report_file_name <- paste0("exp",
                                   parms$exp, "_",
                                   "simOut_",
                                   parms$start_time, 
                                   ".txt")
  parms$results_file_name <- paste0("exp",
                                    parms$exp, "_",
                                    "simOut_", 
                                    parms$start_time,
                                    ".rds")
  parms$description <- c("In each repetition, 1 dataset is created for each condition.
          Imputation methods are used on that condition-specific dataset.
          Results are therefore given per dataset in condition")
  
# Storing prefrences ------------------------------------------------------
  # Parameters to store
  parms$m1_par <- c("rel", "trust.s", "trust.pr") # for model 1
  parms$m2_par <- c("NatAt", "rel", "polAction_r")# for model 2
  
  # Alternative definition
  parms$m1_par <- 1:13
  parms$m2_par <- 1:14
  
  # Needs to match the location and name of the output list
  parms$store <- c(cond         = TRUE,
                   dat_full     = FALSE,
                   dat_miss     = FALSE,
                   m1_EST       = TRUE,
                   m1_CI        = TRUE,
                   m2_EST       = TRUE,
                   m2_CI        = TRUE,
                   fmi          = FALSE,
                   miss_descrps = TRUE,
                   run_time_min = TRUE,
                   imp_values   = FALSE)
  
# Conditions --------------------------------------------------------------
  
  # Fixed random factor
  parms$pm <- c(.1, .2) # make this a random factor
  
  # Experimental factors
  n     <- c(1e3, 3e2) # number of observations
  
  # Dataframe of conditions  
  conds <- expand.grid(n = n)
  
  # Add ridge specification for each condition
  conds <- cbind(conds, ridge = c(1e-2,
                                  1e-07))
  