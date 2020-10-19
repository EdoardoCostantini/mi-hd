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

mod2 <- exp4_fit_mod2(list(GS = dt.f))$GS
cbind(est = coefficients(mod2), confint(mod2))

round(summary(mod2)$adj.r.squared*100, 1)
round(summary(mod2)$coefficients, 3)[, c(1:2, 4)]

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
  