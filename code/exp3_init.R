### Title:    Initialization scirpt, experiment 2
### Project:  Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-07-01

# Fixed Parameters --------------------------------------------------------

  parms <- list()
  parms$exp <- 3 # which experiment is been run

  # Itereations, repetitions, etc
  parms$dt_rep     <- 10# 500 replications for averaging results
  parms$chains     <- 1 # 1   number of parallel chains for convergence check
  parms$iters      <- 5 # 75  total iterations
  parms$burnin_imp <- 0 # 50  how many imputation iterations discarded
  parms$ndt        <- 5 # 10  number of imputed datasets to pool esitmaes from
  parms$thin       <- (parms$iters - parms$burnin_imp)/parms$ndt
  parms$pos_dt  <- (parms$burnin_imp+1):parms$iters # candidates
  parms$keep_dt <- parms$pos_dt[seq(1, 
                                    length(parms$pos_dt), 
                                    parms$thin)] # keep 1 dataset every thin
  
  # For blasso
  parms$chains_bl     <- 1 # 1 
  parms$iters_bl      <- 5 # 2e3  total iterations
  parms$burnin_imp_bl <- 0 # 1950 discarded iterations
  parms$thin_bl       <- (parms$iters_bl - parms$burnin_imp_bl)/parms$ndt
  parms$pos_dt_bl     <- (parms$burnin_imp_bl+1):parms$iters_bl # candidate
  parms$keep_dt_bl    <- parms$pos_dt_bl[seq(1, 
                                             length(parms$pos_dt_bl), 
                                             parms$thin_bl)]
  
  # For mice-like algorithms
  parms$mice_iters <- 5 #  20
  parms$mice_ndt   <- parms$ndt # mice keeps data in a different way
  
# Data gen ----------------------------------------------------------------

  parms$n       <- 200 # number of cases for data generation
  
  # "True" values
  parms$item_mean <- 5
  parms$item_var  <- 5
  
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
  parms$z_m_id  <- c(paste0("z", parms$blck1[1:3]))
  parms$zm_n <- length(parms$z_m_id)
  
  # y gen / imporntant predictors
  parms$yMod_cov <- parms$blck1[1:3]
  parms$yMod_int <- parms$yMod_cov[1:2]
  # parms$formula <- paste0("z", parms$z_m_id, collapse = ", ")
  # parms$formula <- c("y ~ z1 + z2 + z3 + z4 + z5")
  parms$b_main <- rep(1, length(parms$yMod_cov))
  parms$b_int <- 1
  
# Models ------------------------------------------------------------------
  parms$frm <- paste0("y ~ -1 + ",
                      paste0("z", parms$yMod_cov, collapse = " + "))
  parms$frm_int <- paste0("y ~ -1 + ",
                          paste0("z", parms$yMod_cov, collapse = " + "),
                          " + ",
                          paste0("z", parms$yMod_int, collapse = " : ")
                          # paste0("z", parms$yMod_int, collapse = "")
  )
  parms$alphaCI <- .95 # confidence level for parameters CI

# Imputation methods ------------------------------------------------------
  parms$meth_sel <- data.frame(DURR_la = TRUE,
                               DURR_el = FALSE,
                               IURR_la = TRUE,
                               IURR_el = FALSE,
                               blasso  = TRUE,
                               bridge  = TRUE,
                               MI_PCA  = TRUE,
                               MI_CART = TRUE,
                               MI_RF   = TRUE,
                               MI_OP   = TRUE,
                               missFor = TRUE,
                               GS      = TRUE,
                               CC      = TRUE
  )
  parms$methods <- names(parms$meth_sel)[which(parms$meth_sel==TRUE)]
    # (GS, CC always last, alwyas present)
  
  # Location
  parms$missType <- c("high", "low", "tails")[3]
  
  # Response Model (rm)
  # parms$rm_x <- parms$blck1[!parms$blck1 %in% parms$yMod_cov][1:2]
  
  parms$rm_x <- c("y", 
                  paste0("z", 
                         (tail(parms$yMod_cov, 1)+1):
                           (tail(parms$yMod_cov, 1)+2)))
  parms$rm_x_int <- c(parms$rm_x, "int_term")
  
  # weighting the importance of predictors: all the same
  parms$auxWts <- rep(1, length(parms$rm_x))
  parms$auxWts_int <- rep(1, length(parms$rm_x_int))
  
  # IURR
  parms$k_IURR  <- 0 # k value to bias coef sampling covariance matrix
                     # procedure to solve possible issues of singularity

  # PCA
  parms$PCA_inter    <- FALSE # whether you want two way variables interactions
  parms$PCA_poly     <- FALSE # whether you want poly terms 
  parms$PCA_mincor   <- .3 # mincor for qucikpred for single imputation auxiliary vars
  parms$PCA_maxpw    <- 2L # polynomials order
  parms$PCA_pcthresh <- .5 # proportion of vairance for selection of PCs
  parms$formula <- paste0("y", sep = ", ",
                          paste0("z", parms$yMod_cov, collapse = ", "), sep = ", ",
                          paste0(paste0("z", parms$yMod_int), collapse = ""))
  
  
  # Random Forest
  parms$rfntree <- 10
  
  # MICE true
  parms$S_all <- c("y", 
                   paste0("z", parms$blck1),
                   paste0("z", parms$blck2),
                   paste0(paste0("z", parms$yMod_int), collapse = ""))

# Simulation desing -------------------------------------------------------
  # Replicability
  parms$seed     <- 20200512
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
  # Needs to match the location and name of the output list
    
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
  
# Conditions --------------------------------------------------------------

  # Main factors
  int_sub <- c(FALSE, TRUE)
  int_rm  <- c(FALSE, TRUE)
  int_da  <- c(FALSE, TRUE)
  
  # Specifications
  pm    <- c(.3)
  p     <- c(25, 500) # c(50, 500) # number of variables
  r2    <- c(.65)
  
  # Condition dependent Imputation model parameters
  ridge <- 1e-5 # for bridge (needs crossvalidation)
  
  # Dataframe of conditions  
  conds <- expand.grid(pm, ridge, r2, int_sub, int_rm, int_da, p)
  
    colnames(conds) <- c("pm", "ridge", "r2", "int_sub", "int_rm", "int_da", "p")
  conds <- conds[!(conds$p == p[2] & conds$int_da == TRUE), ]
