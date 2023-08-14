# Project:   imputeHD-comp
# Objective: Initialization script for c
# Author:    Edoardo Costantini
# Created:   2022-02-25
# Modified:  2023-04-24
# Notes: 

# Fixed Parameters --------------------------------------------------------

  parms <- list()
  parms$exp <- "1_2" # which experiment is been run

# Iterations, repetitions, etc
  parms$dt_rep     <- 500 # replications for averaging results (500 goal)
  parms$chains     <- 1 # number of parallel chains for convergence check
  parms$iters      <- 70 # 70
  parms$burnin_imp <- 50 # 50 how many imputation iterations should be discarded
  parms$ndt        <- 10 # 10 number of imputed datasets to pool estimates from (10)
  parms$thin       <- (parms$iters - parms$burnin_imp)/parms$ndt
    # every how many iterations should you keep the imputation for a dataset
    # Example: of 20 iterations, I burn the first 10 I need for convergence
    # then, I keep 1 set of imputations every "thin"
    # number of iterations for the MI procedure (20)
  parms$pos_dt  <- (parms$burnin_imp+1):parms$iters # candidate datasets (after convergence)
  parms$keep_dt <- parms$pos_dt[seq(1, length(parms$pos_dt), parms$thin)] # keep 1 dataset every thin

  # For blasso
  parms$chains_bl     <- 1 # number of parallel chains for convergence check
  parms$iters_bl      <- 70 # 70
  parms$burnin_imp_bl <- 50 # 50 how many imputation iterations should be discarded
  parms$thin_bl       <- (parms$iters_bl - parms$burnin_imp_bl)/parms$ndt
  parms$pos_dt_bl     <- (parms$burnin_imp_bl+1):parms$iters_bl # candidate datasets
  parms$keep_dt_bl    <- parms$pos_dt_bl[seq(1, 
                                             length(parms$pos_dt_bl), 
                                             parms$thin_bl)]
  
  # For mice-like algorithms
  parms$mice_iters <- 50 # 50
  parms$mice_ndt   <- parms$ndt # 10 # number of imputed datasets to pool estimates from (10)
  
# Data gen ----------------------------------------------------------------

# Sample size for data
  parms$n       <- 200

# "True" values of variables
  parms$item_mean <- 5
  parms$item_var  <- 5
  
# Variable blocks
  parms$blck1 <- 1:5 # highly correlated
  parms$blck2 <- (tail(parms$blck1,1)+1):10 # correlated

  parms$blck1_r <- .6 # correlation for highly correlated variables
  parms$blck2_r <- .3 # correlation for correlated variables

# Z gen (fully observed covariates)
  parms$Z_o_mu  <- 0 # mean for gen of fully observed variables

# z_m_id gen (variables that will have missingness)
  parms$z_m_id <- paste0("z", c((parms$blck1[1]):3, (parms$blck2[1]):8))
  parms$zm_n <- length(parms$z_m_id)
  parms$S_all <- list(
    q1 = (paste0("z", c(parms$blck1, parms$blck2))), # variables in block 1 and 2
    q2 = (c(4:13, 50:59))
  )[[1]]

# y gen / important predictors
  parms$formula <- paste0(parms$z_m_id, collapse = ", ")
  parms$lm_model <- paste0(parms$z_m_id[1], 
                           " ~ - 1 + ", 
                           paste0(parms$z_m_id[-1], collapse = " + "))

# Response Model (rm)
  parms$missType <- c("high", "low")[2]
  parms$auxWts <- c(1, 1, 1, 1) # weighting the importance of predictors: all the same
  parms$rm_x <- matrix(
    c(
      4, 5, 9, 10,
      4, 5, 9, 10,
      4, 5, 9, 10,
      4, 5, 9, 10,
      4, 5, 9, 10,
      4, 5, 9, 10
    ),
    ncol = 4, nrow = 6,
    byrow = TRUE
  )

# Models ------------------------------------------------------------------
  # source("./gen_lavaan_model.R") # generate txt file for lavaan model
  ## NOTE: this is not to be kept running in the windows simulation because if so
  ##       it creates problems. Run it just once to have it in the correct folder.
  parms$lav_model <- paste(readLines("../txt/lavaan_model_sat.txt"), collapse="\n")
  
  # Methods
  parms$meth_sel <- data.frame(
    DURR_la = FALSE,
    IURR_la = FALSE,
    blasso = FALSE,
    bridge = FALSE,
    MI_PCA = FALSE,
    MI_CART = FALSE,
    MI_RF = FALSE,
    stepFor = TRUE,
    MI_qp = FALSE,
    MI_am = FALSE,
    MI_OP = FALSE,
    missFor = TRUE,
    mean = TRUE,
    CC = TRUE,
    GS = TRUE
  )
  parms$methods <- names(parms$meth_sel)[which(parms$meth_sel==TRUE)]
    # (GS, CC always last, always present)

# Analysis related
  parms$alphaCI <- .95 # confidence level for parameters CI

# IURR related
  parms$k_IURR  <- 0 # k value to bias coef sampling covariance matrix in IURR
                     # procedure to solve possible issues of singularity

# PCA method related
  parms$PCA_pcthresh <- .5 # proportion of variance for selection of PCs

# Random forest related
  parms$rfntree <- 2  # default value as no indication is given in paper
  parms$missFor_maxiter <- 2 # maxiter = 20
  parms$missFor_ntree <- 2  # ntree = 100

# Detect computer OS
parms$OS <- switch(Sys.info()[["sysname"]],
  Windows = {
    print("Windows")
  },
  Linux = {
    print("Linux")
  },
  Darwin = {
    print("Mac")
  }
)

# IVE specifics
parms$IVEloc <- c(
  Windows = "C:\\Program Files\\Srclib\\R", # or other location https://www.src.isr.umich.edu/software/iveware/iveware-documentation/installation-guide/
  Mac = "/Library/Srclib/R"
)[parms$OS]

# Replicability related
  parms$seed     <- 20230428
  parms$nStreams <- 1000

# Storing preferences ------------------------------------------------------
  # Needs to match the location and name of the output list
 
  parms$store <-  c(cond         = TRUE,
                    dat_full     = FALSE,
                    dat_miss     = FALSE,
                    sem_EST      = TRUE,
                    sem_CI       = TRUE,
                    lm_EST       = TRUE,
                    lm_CI        = TRUE,
                    fmi          = FALSE,
                    miss_descrps = FALSE,
                    run_time_min = TRUE,
                    imp_values   = FALSE)
  
# Output and Progress report related
  parms$outDir <- "../output/"
  parms$start_time <- format(Sys.time(), "%Y%m%d_%H%M")
  parms$report_file_name <- paste0(
    "exp",
    parms$exp, "_",
    "simOut_",
    parms$start_time,
    ".txt"
  )
  parms$results_file_name <- paste0(
    "exp",
    parms$exp, "_",
    "simOut_",
    parms$start_time,
    ".rds"
  )
  parms$description <- c("In each repetition, 1 dataset is created for each condition.
        Imputation methods are used on that condition-specific dataset.
        Results are therefore given per dataset in condition")

# Conditions --------------------------------------------------------------

  # Collinearity study: experimental factor levels
  pm <- .3
  p   <- c(50, 500)
  latent <- FALSE
  collinearity <- c(NA, .6, .8, .9)
  ridge <- c(1e-1, 1e-07, rep(1e-8, 3), 1e-7, rep(1e-6, 2))
  minR2 <- 1e-3

  # Create experimental conditions
  conds <- expand.grid(
    p = p,
    latent = latent,
    pm = pm,
    collinearity = collinearity
  )

  # Sort based on p
  conds <- dplyr::arrange(conds, p)

  # Bridge special parameters per condition
  conds$ridge <- ridge

  # IVEware special parameters per condition
  conds$minR2 <- minR2