### Title:    Initialization scirpt parameter for experiment 3 (simple latent)
### Project:  Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-07-22

# Fixed Parameters --------------------------------------------------------

  parms <- list()
  parms$exp <- 2 # second experiment - latent structure

# Itereations, repetitions, etc
  parms$dt_rep     <- 5 # 500 replications for averaging results (200 goal)
  parms$chains     <- 1 # 1   number of parallel chains for convergence check
  parms$iters      <- 5 # 75
  parms$burnin_imp <- 0 # 50  how many imputation iterations should be discarded
  parms$ndt        <- 5 # 10  number of imputed datasets to pool esitmaes from (10)
  parms$thin       <- (parms$iters - parms$burnin_imp)/parms$ndt
    # every how many iterations should you keep the imputation for a dataset
    # Example: of 20 iterations, I burn the first 10 I need for convergence
    # then, I keep 1 set of imputations every "thin"
    # number of iterations for the MI procedure (20)
  parms$pos_dt  <- (parms$burnin_imp+1):parms$iters # candidate datasets (after convergence)
  parms$keep_dt <- parms$pos_dt[seq(1, 
                                    length(parms$pos_dt), 
                                    parms$thin)] # keep 1 dataset every thin

  # For blasso
  parms$chains_bl     <- 1 # 1 number of parallel chains for convergence check
  parms$iters_bl      <- 5 # 2e3
  parms$burnin_imp_bl <- 0 # 1950 how many imputation iterations should be discarded
  parms$thin_bl       <- (parms$iters_bl - parms$burnin_imp_bl)/parms$ndt
  parms$pos_dt_bl     <- (parms$burnin_imp_bl+1):parms$iters_bl # candidate datasets
  parms$keep_dt_bl    <- parms$pos_dt_bl[seq(1, 
                                             length(parms$pos_dt_bl), 
                                             parms$thin_bl)]
  
  # For mice-like algorithms
  parms$mice_iters <- 5 #  20
  parms$mice_ndt   <- parms$ndt # 10 # number of imputed datasets to pool esitmaes from (10)
  
# Data gen ----------------------------------------------------------------
  parms$n_obs <- 200 # number of cases for data generation
  parms$n <- 200 # number of cases for data generation
  
# "True" values
  parms$lv_mean   <- 0
  parms$lv_var    <- 1
  parms$item_mean <- 5 # Mean not zero so that bias in percent makes more sense
  parms$item_var  <- 1
  
# Latent Structure
  parms$n_it    <- 5 # number of measured items
  
# Latent Variable blocks
  parms$blck1 <- 1:4 # highly correlated latent variables
  parms$blck2 <- (tail(parms$blck1, 1)+1):8 # mid correlated latent variables

  parms$blck1_r <- .6 # correlation for highly correlated variables
  parms$blck2_r <- .3 # correlation for correlated variables

# Z gen (fully observed covariates)

# z_m gen (measured items that will have missingness)
  # parms$z_m_id  <- 1:10
  parms$z_m_id <- paste0("z", 1:10)
  parms$zm_n <- length(parms$z_m_id)
  parms$S_all   <- paste0("z", ( 1:(parms$n_it * tail(parms$blck1, 1)) ))
    # all measured items (5) for the first 4 lv. These include:
    # - latent variables of items with missing values (1 and 2)
    # - latent variables involved in response model (3 and 4)
    # - all of their items
  
# y gen / imporntant predictors
  # parms$formula <- paste0("z", parms$z_m_id, collapse = ", ")
  parms$formula <- paste0(parms$z_m_id, collapse = ", ")
  # variables that are to be imputed
  # parms$lm_model <- paste0("z", parms$z_m_id) # not in formula version
  parms$lm_model <- parms$z_m_id
  
# Response Model
  parms$missType <- c("high", "low")[2]
  
  # Latent variables 3 and 4 influence missingness
  parms$rm_x <- (2*parms$n_it+1):(4*parms$n_it)
  
  # weighting the importance of predictors: all the same
  parms$auxWts <- rep(1, length(parms$rm_x))
  
# Models ------------------------------------------------------------------
  parms$lav_model <- NULL #paste(readLines("../txt/lavaan_model_sat.txt"), 
                          #collapse="\n")
  parms$sc_n <- 3 # how many "Scores" in the sat model for SCore data
  
# Generic
  parms$meth_sel <- data.frame(DURR_la    = TRUE,
                               IURR_la    = TRUE,
                               blasso     = TRUE,
                               bridge     = TRUE,
                               MI_PCA     = TRUE,
                               MI_CART    = TRUE,
                               MI_RF      = TRUE,
                               MI_OP      = TRUE,
                               missFor    = TRUE,
                               mean       = TRUE,
                               CC         = TRUE,
                               GS         = TRUE)
  parms$methods <- names(parms$meth_sel)[which(parms$meth_sel==TRUE)]
    # (GS, CC always last, alwyas present)

  parms$alphaCI <- .95 # confidence level for parameters CI
  parms$k_IURR  <- 0 # k value to bias coef sampling covariance matrix in IURR
                     # procedure to solve possible issues of singularity

# PCA method related
  parms$PCA_inter    <- FALSE # whether you want two way variables interactions
  parms$PCA_poly     <- FALSE # whether you want poly terms 
  parms$PCA_mincor   <- .3 # mincor for qucikpred for single imputation auxiliary vars
  parms$PCA_maxpw    <- 2L # polynomials order
  parms$PCA_pcthresh <- .5 # proportion of vairance for selection of PCs

# Random forest related
  parms$rfntree <- 10  # default value as no indication is given in paper
  parms$missFor_maxiter <- 20 # maxiter = 20
  parms$missFor_ntree <- 100  # ntree = 100
  
# Bayesina Ridge imputation related
  parms$ridge <- 1e-5
  
# Replicability related
  parms$seed     <- 20200512 #20200309
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
  parms$store <- c(cond     = TRUE,
                   dat_full = FALSE,
                   dat_miss = FALSE,
                   sem_EST  = TRUE,
                   sem_CI   = TRUE,
                   CFA_EST  = TRUE,
                   CFA_CI   = TRUE,
                   semS_EST = TRUE,
                   semS_CI  = TRUE,
                   lm_EST   = TRUE,
                   lm_CI    = TRUE,
                   fmi      = TRUE,
                   miss_des = FALSE,
                   time     = TRUE,
                   imps     = FALSE)
  
# Conditions --------------------------------------------------------------

  # Define Experimental Factor Values
  lv    <- c(10, 100)       # number of latent variables
  pm    <- c(.1, .3)        # proportion of missings level
  fl    <- c("high", "low") # factor loadings level
  ridge <- rep(c(1e-1, 1e-7), 4) # 1 valude found w/ corssvalidation
  
  # Make Conditions
  conds <- expand.grid(lv, pm, fl,
                       stringsAsFactors = FALSE)
  
  # Append Ridge Parameter
  conds <- cbind(conds, ridge)
  
  # Select Conditions for run
  conds <- conds[1:(nrow(conds) - 0), ]
  
  # Give Meaningful names to Columns
  colnames(conds) <- c("lv", "pm", "fl", "ridge")

  # Print a latex table    
  xtable::xtable(conds[, c("fl", "pm", "lv")])
  