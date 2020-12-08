### Title:    Checking parts of experiment 4 work as expected
### Project:  Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-10-06

# Substantive models ------------------------------------------------------

rm(list=ls())
source("./init_general.R")
source("./exp4_init.R")

EVS_dt <- readRDS("../data/exp4_EVS2017_full.rds")

dt.o <- EVS_dt$orig
dt.f <- EVS_dt$full

mod1 <- exp4_fit_mod1(list(GS = dt.f))$GS
cbind(est = coefficients(mod1), confint(mod1))

round(summary(mod1)$adj.r.squared*100, 1)
round(summary(mod1)$coefficients, 3)[, c(1:2, 4)]
nrow(summary(mod1)$coefficients)

mod2 <- exp4_fit_mod2(list(GS = dt.f))$GS
cbind(est = coefficients(mod2), confint(mod2))

round(summary(mod2)$adj.r.squared*100, 1)
round(summary(mod2)$coefficients, 3)[, c(1:2, 4)]
nrow(summary(mod2)$coefficients)

# CIR Full data -----------------------------------------------------------

  rm(list=ls())
  source("./init_general.R")
  source("./exp4_init.R")

# Select 1 condition
  cond <- conds[2, ]
  
# get dataset
  set.seed(20200805)
  data_source <- readRDS("../data/exp4_EVS2017_full.rds")$full
  reps <- 1e3
  m1_CI <- matrix(0, ncol = 2, nrow = 13)
  m2_CI <- matrix(0, ncol = 2, nrow = 14)

# Store results
  mpars <- cov_o <- list(m1 = matrix(NA, nrow = reps, ncol = 13),
                         m2 = matrix(NA, nrow = reps, ncol = 14))
  CI <- NULL
# Perform check
for (r in 1:reps) {
  # Gen one fully-obs data
  Xy <- data_source[sample(1:nrow(data_source),
                           cond$n,
                           replace = TRUE), ]
  # Fit model 1 and 2
  si_data <- list(GS = Xy)
  
  ## Fit --------------------------------------------------------------------- ##
  # model 1
  m1_sn <- exp4_fit_mod1(si_data)
    mpars[[1]][r,]  <- lapply(m1_sn, coef)$GS
  m2_sn <- exp4_fit_mod2(si_data)
    mpars[[2]][r,]  <- lapply(m2_sn, coef)$GS
  
  # Get Confidence intervals
  CI[[r]] <- lapply(c(m1 = m1_sn, m2 = m2_sn), confint)
}
  
  # Obtain CIR results
  mpar_true <- lapply(mpars, colMeans)
  for (r in 1:reps) {
    cov_o[[1]][r, ] <- CI[[r]]$m1.GS[, 1] < mpar_true$m1 & mpar_true$m1 < CI[[r]]$m1.GS[, 2]
    cov_o[[2]][r, ] <- CI[[r]]$m2.GS[, 1] < mpar_true$m2 & mpar_true$m2 < CI[[r]]$m2.GS[, 2]
  }
  colnames(cov_o[[1]]) <- names(lapply(m1_sn, coef)$GS)
  colnames(cov_o[[2]]) <- names(lapply(m2_sn, coef)$GS)
  
  lapply(cov_o, colMeans)$m1["rel"]

# Identify Variables for MAR ----------------------------------------------
# Criteria:
# a) correlated with the target of missing variables
# b) not target of missing
# c) not all part of the substantive models
# d) they make sense as response mechanisms

  rm(list=ls())
  source("./init_general.R")
  source("./exp4_init.R")
  cond <- conds[1, ]
  
  # get dataset
  set.seed(20200805)
  data_source <- readRDS("../data/exp4_EVS2017_full.rds")$full

  # Consider correlation between target variables and all predictors
  mm <- model.matrix(~., data_source)
  cormm <- cor(mm)[, parms$z_m_id]
  
  # Which of them are higher than .1
  store <- NULL
  for (i in 1:6) {
    store[[i]] <- names(which(abs(cormm[, i]) >= .1))
  }
  
  # Which possible preidctors are correlated .1 or + with all targets
  v_presence <- NULL
  for (v in 1:ncol(mm)) {
    pat <- paste0("\\<", colnames(mm)[v], "\\>")
    v_presence[[v]] <- grep(pat, store)  
    
  }
  
  # get the names
  possible <- colnames(mm)[which(sapply(v_presence, length) == 6)]
  cormm[possible, ]
  
  # v35: trust in people you just met is the only one (other than edu)
  # that makes sense as a response predictor: if you don't trust the 
  # interviewer you just met, you will maybe not disclose information
  
# Effects of missingness imposition on estiamtes --------------------------
# a) Saturated model on raw data (only items with missing values)
# b) CFA model on raw data (latent variables target of miss and 
#    1 latent variable influencing)
# c) Saturated Model on Scores obtained from the variables (mean of items)
# d) LM on Scores
  
  rm(list=ls())
  source("./init_general.R")
  source("./exp4_init.R")
  
  # Tweak paramters
  cond <- conds[1, ]
  parms$m1_par <- c("rel", "trust.s", "trust.pr")
  parms$m2_par <- c("rel", "NatAt", "polAction_r")
  
  # Create storing objects
  reps <- 1e3
  output <- vector("list", reps)
  # m1_est <- m2_est <- data.frame(MEAN = NA, CC = NA, GS = NA)
  # m1_est <- list()
  # m2_est <- data.frame(MEAN = NA, CC = NA, GS = NA)
  # m1_ci <- m2_ci <- list()
  set.seed(20200805)
  data_source <- readRDS("../data/exp4_EVS2017_full.rds")$full
  
  # True Values
  b.true <- list(m1 = lm_EST(exp4_fit_mod1(list(GS = data_source)))[parms$m1_par,],
                 m2 = lm_EST(exp4_fit_mod2(list(GS = data_source)))[parms$m2_par,])
  
  # Perform analysis
  pb <- txtProgressBar(min = 0, max = reps, style = 3)
  for (r in 1:reps) {
    output[[r]] <- list(cond1 = list(), cond2 = list())
    # Gen one fully-obs data
    Xy <- data_source[sample(1:nrow(data_source), 
                             cond$n,
                             replace = TRUE), ]
    # Impose missing values
    Xy_mis <- imposeMiss_evs(Xy, parms, cond)
    # Missing data
    miss_descrps <- colMeans(is.na(Xy_mis)[, parms$z_m_id])
  
  # Mean imputation  
    Xy_mean <- Xy_mis
    for (j in 1:parms$zm_n) {
      ry <- !is.na(Xy_mean[, parms$z_m_id[j]])
      mu_zj <- mean(Xy_mean[, parms$z_m_id[j]], na.rm = TRUE)
      Xy_mean[!ry, parms$z_m_id[j]] <- mu_zj
    }
    
  # Fit models to single datasets
    si_data <- list(mean    = Xy_mean,
                    CC      = na.omit(Xy_mis),
                    GS      = Xy)
  # Fit 
    # Model 1
    m1_sn <- exp4_fit_mod1(si_data)
    # Model 2
    m2_sn <- exp4_fit_mod2(si_data)
    
  # Extraction indeces
    m1_EST_indx <- parms$m1_par
    m1_CI_indx <- c(do.call(rbind, lapply(m1_EST_indx, grep,
                                          rownames(lm_CI(m1_sn)))))
    m2_EST_indx <- parms$m2_par
    m2_CI_indx <- c(do.call(rbind, lapply(m2_EST_indx, grep,
                                          rownames(lm_CI(m2_sn)))))
    
  # Store Output
    output[[r]][[1]] <- list(# Model 1
      m1_EST       = lm_EST(m1_sn)[m1_EST_indx, ],
      m1_CI        = lm_CI(m1_sn)[m1_CI_indx, ],
      # Model 2
      m2_EST       = lm_EST(m2_sn)[m2_EST_indx, ],
      m2_CI        = lm_CI(m2_sn)[m2_CI_indx, ])

    # Monitor progress
    setTxtProgressBar(pb, r)
  }
  
  output
  output$parms$methods <- colnames(output[[1]]$cond1$m1_EST)
  output$parms$dt_rep  <- reps
  output$parms$exp     <- 4
  
  # Obtain results
  m1_res <- res_sum(output, 
                    model = "m1", 
                    condition = 1)
  cbind(bt = round(b.true$m1, 3), m1_res$bias_per)
  m1_res$ci_cov
  
  m2_res <- res_sum(output, 
                    model = "m2", 
                    condition = 1)
  cbind(bt = round(b.true$m2, 3), m2_res$bias_per)
  m2_res$ci_cov
  
# Optimal Model imputation ------------------------------------------------

  rm(list=ls())
  source("./init_general.R")
  source("./exp4_init.R")
  
  ## Define inputs --------------------------------------------------------- ##
  
  set.seed(1234)
  cond <- conds[2, ]
  data_source <- readRDS("../data/exp4_EVS2017_full.rds")$full
  
  ## Data ------------------------------------------------------------------ ##
  # Gen one fully-obs data
  Xy <- data_source[sample(1:nrow(data_source), 
                           cond$n,
                           replace = TRUE), ]
  
  # Impose missing values
  Xy_mis <- imposeMiss_evs(Xy, parms, cond)
  
  # Missing data
  miss_descrps <- colMeans(is.na(Xy_mis)[, parms$z_m_id])
  
  ## Impute ---------------------------------------------------------------- ##
  
  # Impute with "Optimal model"
  imp_MICE_OP <- impute_MICE_OP(Z = Xy_mis,
                                O = data.frame(!is.na(Xy_mis)),
                                cond = cond,
                                perform = parms$meth_sel$MI_OP,
                                parms = parms)
  plot(imp_MICE_OP$mids)
  
  # Inside this happens:
  Z = Xy_mis
  O = data.frame(!is.na(Z))
  Z[, parms$z_m_id]
  cond = cond
  perform = parms$meth_sel$MI_OP
  parms = parms
  
  # Define predictor matrix for MI TRUE with best active set

  # Try with pred matrix
  predMat <- quickpred(Z,
                       mincor = .3,
                       include = parms$S_all)
  rowSums(predMat[p_imp_id, ])

  # Impute
  imp_MITR_mids <- mice::mice(Z, 
                              predictorMatrix = predMat,
                              m = parms$mice_ndt,
                              maxit = parms$mice_iters,
                              ridge = 1e-5,
                              method = "pmm")
  imp_MITR_mids$loggedEvents[1:5]
  plot(imp_MITR_mids)
  

# Include all predictors --------------------------------------------------
# The proble of including all predictors
  rm(list=ls())
  source("./init_general.R")
  source("./exp4_init.R")
  
  ## Define inputs --------------------------------------------------------- ##
  
  set.seed(1234)
  cond <- conds[1, ]
  data_source <- readRDS("../data/exp4_EVS2017_full.rds")$full
  
  ## Data ------------------------------------------------------------------ ##
  # Gen one fully-obs data
  Xy <- data_source[sample(1:nrow(data_source), 
                           cond$n,
                           replace = TRUE), ]
  
  # Impose missing values
  Xy_mis <- imposeMiss_evs(Xy, parms, cond)
  
  # Missing data
  miss_descrps <- colMeans(is.na(Xy_mis)[, parms$z_m_id])
  
  ## Impute ---------------------------------------------------------------- ##
  
  # Impute with "Optimal model"
  imp_MICE_OP <- impute_MICE_OP(Z = Xy_mis,
                                O = data.frame(!is.na(Xy_mis)),
                                cond = cond,
                                perform = parms$meth_sel$MI_OP,
                                parms = parms)
  plot(imp_MICE_OP$mids)
  
  # Inside this happens:
  Z = Xy_mis
  O = data.frame(!is.na(Z))
  Z[, parms$z_m_id]
  cond = cond
  perform = parms$meth_sel$MI_OP
  parms = parms
  
  # Define predictor matrix for MI TRUE with best active set
  
  # Include all predictors
  dim(Z)
  predMat <- make.predictorMatrix(Z)
  dim(predMat)

  # Impute
  pdraw <- .norm.draw(y       = target, 
                      ry      = rj,
                      x       = Zx,
                      ls.meth = "ridge", 
                      ridge   = 1)
  estimice
  pdraw <- .norm.draw(y       = target, 
                      ry      = rj,
                      x       = Zx,
                      ls.meth = "non-ridge method", 
                      ridge   = ridge_p)
  
  imp_MITR_mids <- mice::mice(Z, 
                              predictorMatrix = predMat,
                              m = parms$mice_ndt,
                              maxit = parms$mice_iters,
                              ridge = 1,
                              method = "pmm")
  imp_MITR_mids$loggedEvents[1:5]
  plot(imp_MITR_mids)

# Manually
  
  p  <- ncol(Z) # number of variables [indexed with J]
  
  p_imp    <- sum(colMeans(O) < 1)
  p_imp_id <- names(which(colMeans(O) < 1))
  nr       <- colSums(!O[, colMeans(O) < 1])
  
  # To store imputed values and check convergence
  imp_bridge_val <- vector("list", parms$chains)
  names(imp_bridge_val) <- seq(1:parms$chains)
  
  # Time performance
  start.time <- Sys.time()  

  # To store multiply imputed datasets (in from the last chain)
  imp_bridge_dat <- vector("list", parms$iters)
  names(imp_bridge_dat) <- seq(1:parms$iters)
  
  Zm <- init_dt_i(Z, missing_type(Z)) # initialize data for each chain
  imp_bridge_dat$`1` <- Zm
  # Empty storing objects for MCMC samples
  
  # Imputed scores
  Imp.out <- lapply(p_imp_id, function(x) {
    matrix(data = NA, nrow = parms$iters, ncol = nr[x],
           dimnames = list(NULL, rownames(Zm[!O[, x],]) ))
  })
  for (i in 1:p_imp) Imp.out[[i]][1, ] <- Zm[!O[, p_imp_id[i]], 
                                             p_imp_id[i]]
  
  for (m in 2:parms$iters) {
    # Loop across variables (cycle)
    for (j in 1:p_imp) {
      j <- 1
      J <- which(colnames(Zm) %in% p_imp_id[j])
      rj <- !is.na(Z[, J])
      zj_obs <- Zm[rj, J]
      zj_mis <- Zm[!rj, J]
      
      # PREP data for post draw
      target   <- Zm[, p_imp_id[j]]
      Zx       <- model.matrix(~ ., Zm[, -J])[, -1]
      
      # Find Constants
      const <- names(which(apply(Zx, 2, var) == 0))
      Zx    <- Zx[, !colnames(Zx) %in% const]
      
      # Find Dummies that have 99% of objservations in 1 category
      tabular <- apply(Zx, 2, table)
      
      # Select only dummies
      tabular <- tabular[sapply(tabular, length) == 2]
      
      # Vector of dummy names to discard
      dum.disc <- lapply(tabular, function(x) {
        x[1] / sum(x) > .95 | x[1] / sum(x) < 1-.95
      })
      dum.disc <- names(which(dum.disc == TRUE))
      Zx    <- Zx[, !colnames(Zx) %in% dum.disc]
      
      # Find collinear variables
      coll.vars <- find.collinear(Zx)
      Zx  <- Zx[, !colnames(Zx) %in% coll.vars]
      
      # Obtain Post Draw
      pdraw <- .norm.draw(y       = target, 
                          ry      = rj,
                          x       = Zx,
                          ls.meth = "ridge", 
                          ridge   = 0)
      
      p <- estimice(Zx[rj, ], target[rj])
      x = Zx[rj, 1:50]
      y = target[rj]
      ridge = 0
      xtx <- crossprod(x)
      pen <- ridge * diag(xtx)
      if (length(pen) == 1) 
        pen <- matrix(pen)
      v <- solve(xtx + diag(pen))
      c <- t(y) %*% x %*% v
      r <- y - x %*% t(c)
      return(list(c = t(c), r = r, v = v, df = df, ls.meth = ls.meth))
      
      # Fix Z_mis to drop columns that are not relevant anymore
      Z_mis <- Zm[!rj, -J]
      Z_mis <- model.matrix(~ ., Z_mis)[, -1]
      Z_mis <- Z_mis[, !colnames(Z_mis) %in% c(const, dum.disc, coll.vars)]
      
      # Obtain posterior predictive draws
      pdraw_zj_imp <- Z_mis %*% pdraw$beta + rnorm(sum(!rj)) * pdraw$sigma
      
      # Append imputation (for next iteration)
      Zm[!rj, J] <- pdraw_zj_imp
      
      # Store imputations
      Imp.out[[j]][m, ]  <- pdraw_zj_imp
    }
    # only data from last chain will actually be saved
    imp_bridge_dat[[m]] <- Zm
  }
  imp_bridge_val[[cc]] <- Imp.out
