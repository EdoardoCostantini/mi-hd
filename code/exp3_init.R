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
  parms$iters      <- 2 # 75  total iterations
  parms$burnin_imp <- 0 # 50  how many imputation iterations discarded
  parms$ndt        <- 2 # 10  number of imputed datasets to pool esitmaes from
  parms$thin       <- (parms$iters - parms$burnin_imp)/parms$ndt
  parms$pos_dt  <- (parms$burnin_imp+1):parms$iters # candidates
  parms$keep_dt <- parms$pos_dt[seq(1, 
                                    length(parms$pos_dt), 
                                    parms$thin)] # keep 1 dataset every thin
  
  # For blasso
  parms$chains_bl     <- 1 # 1 
  parms$iters_bl      <- 2 # 300  total iterations
  parms$burnin_imp_bl <- 0 # 250 discarded iterations
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
  parms$z_m_id  <- c(paste0("z", parms$blck1[1:3]))
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
  parms$meth_sel <- data.frame(DURR_la    = TRUE,
                               DURR_SI    = TRUE,
                               IURR_la    = TRUE,
                               IURR_SI    = TRUE,
                               blasso     = TRUE,
                               blasso_SI  = TRUE,
                               bridge     = TRUE,
                               bridge_SI  = TRUE,
                               MI_PCA     = TRUE,
                               MI_CART    = TRUE,
                               MI_CART_SI = TRUE,
                               MI_RF      = TRUE,
                               MI_RF_SI   = TRUE,
                               MI_OP      = TRUE,
                               missFor    = TRUE,
                               GS         = TRUE,
                               CC         = TRUE)
  parms$methods <- names(parms$meth_sel)[which(parms$meth_sel==TRUE)]
    # (GS, CC always last, alwyas present)
  
  # Location
  parms$missType <- c("high", "low", "tails")[3]
  
  # Response Model (rm)
  # parms$rm_x <- c("y", 
  #                 paste0("z", 
  #                        (tail(parms$yMod_cov, 1)+1):
  #                          (tail(parms$yMod_cov, 1)+2)))
  parms$rm_x <- c("y", 
                  paste0("z", 
                         (tail(parms$yMod_cov, 1))))
  parms$rm_x_int <- c(parms$rm_x, "int_term")
  
  # weighting the importance of predictors: all the same
  parms$auxWts <- rep(1, length(parms$rm_x))
  parms$auxWts_int <- rep(1, length(parms$rm_x_int))
  
  # IURR
  parms$k_IURR  <- 0 # k value to bias coef sampling covariance matrix
                     # procedure to solve possible issues of singularity

  # PCA
  parms$SI_iter      <- 10L   # 200L # iterations for single imputation in PCA run

  # parms$formula <- paste0("y", sep = ", ",
  #                         paste0("z", parms$yMod_cov, collapse = ", "), sep = ", ",
  #                         paste0(paste0("z", parms$yMod_int), collapse = ""))
  
  # Random Forest
  parms$rfntree <- 10
  parms$missFor_maxiter <- 20 # maxiter = 20
  parms$missFor_ntree <- 100  # ntree = 100
  
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
                   fmi          = FALSE,
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
  p     <- c(25, 500) # c(25, 500) # number of variables
  r2    <- c(.5)
    # Note that if you are excluding the intercept you are inflating the
    # R squared. Howver, this does not really matter for the inference on
    # the regression coefficients.
  
  # Condition dependent Imputation model parameters
  ridge <- c(1e-08, # values from the exp3_cv_bridge_results.R script
             1e-01, 
             1e-04, 
             1e-01,
             1e-08, 
             1e-08, 
             1e-07, 
             1e-07, 
             1e-07, 
             1e-07, 
             1e-07, 
             1e-07)
  
  # Dataframe of conditions  
  conds <- expand.grid(pm, ridge=NA, r2, int_sub, int_rm, int_da, p)
  
    colnames(conds) <- c("pm", "ridge", "r2", "int_sub", "int_rm", "int_da", "p")
  conds <- conds[!(conds$p == p[2] & conds$int_da == TRUE), ]
  conds$ridge <- ridge
