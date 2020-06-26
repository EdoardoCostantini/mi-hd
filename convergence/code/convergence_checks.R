### Title:    performing pre-sim convergence checks
### Project:  Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-05-19
### Notes:    The goal is to find out how many iterations should be used
###           in the full study for each method.

rm(list=ls())
source("./init.R")
source("./functions.R")

# Location for saving imputations for convergence check
path_loc <- "../convergence/output/"

# Run of all models of interest
parms$meth_sel <- data.frame(DURR_rd = TRUE,
                             DURR_la = TRUE,
                             DURR_el = TRUE,
                             IURR_la = TRUE,
                             IURR_el = TRUE,
                             blasso  = TRUE,
                             MI_PCA  = TRUE,
                             MI_CART = TRUE,
                             MI_CTBB = FALSE,
                             MI_RF   = TRUE,
                             MI_T    = TRUE,
                             missFor = TRUE,
                             GS      = TRUE,
                             CC      = TRUE
)

# Dimensionality conditions to check
# Check  convergence for most difficult condition
# and use decision for all conditions

p   <- c(500) # number of variables
rho <- c(0.5) # autoregressive structure
q   <- c(1) # c(4, 20)  # active set (num of variables are true predictors of y)
latent <- c(FALSE, TRUE)[1]
pm <- c(.3)

cond_convcheck <- expand.grid(p, rho, q, latent, pm)
  colnames(cond_convcheck) <- c("p", "rho", "q", "latent", "pm")
  
cond_convcheck
  
# Get one dataset
  Xy <- genDt_mvn(cond_convcheck[1,], parms)
  Xy_mis <- imposeMiss(Xy, parms, cond_convcheck[1,])
  Xy_mis <- cbind(Xy_mis[parms$z_m_id], Xy_mis[-parms$z_m_id])
  
  O <- !is.na(Xy_mis) # matrix index of observed values
  miss_descrps <- colMeans(!O[, 1:parms$zm_n]) 

# Check Convergences One Method at the time -------------------------------
  
# DURR lasso --------------------------------------------------------------
  
  parms$chains     <- 5 # number of parallel chains for convergence check
  parms$iters      <- 100
  parms$burnin_imp <- 0#10 # how many imputation iterations should be discarded
  parms$ndt        <- 1#10 # number of imputed datasets to pool esitmaes from (10)
  parms$thin       <- (parms$iters - parms$burnin_imp)/parms$ndt
  parms$pos_dt  <- (parms$burnin_imp+1):parms$iters # candidate datasets (after convergence)
  parms$keep_dt <- parms$pos_dt[seq(1, length(parms$pos_dt), parms$thin)] # keep 1 dataset every thin
  imp_DURR_la <- impute_DURR(Z = Xy_mis,
                             O = as.data.frame(O),
                             reg_type = "lasso",
                             perform = parms$meth_sel$DURR_la,
                             parms = parms)
  
  # Store/load results
  saveRDS(object = list(imp_obj = imp_DURR_la,
                        parms = parms,
                        cond = cond_convcheck),
          file = paste0(path_loc, "DURR_cnvg", ".rds") )
  
  imp_blasso <- readRDS(paste0(path_loc, "DURR_cnvg", ".rds"))$imp_obj
  
  # Convergence Plot
  par(mfrow = c(2,3))
  for (v in 1:length(parms$z_m_id)) {
    imps_4plot <- imp_DURR_la$imps[[1]][[v]]
    mean_imp <- rowMeans(imps_4plot)
    plot(1:parms$iters, mean_imp, type = "l",
         main = "Blasso Mean Imputations",
         ylim = c(mean(mean_imp)-4*sd(mean_imp), mean(mean_imp)+4*sd(mean_imp)),
         ylab = paste0("z", v), xlab = "Iteration")
    # And add imputations from other chains
    for (i in 2:(parms$chains)) {
      imps_4plot <- imp_DURR_la$imps[[i]][[v]]
      mean_imp <- rowMeans(imps_4plot)
      lines(1:parms$iters, mean_imp)
    }
  }
  
# Blasso ------------------------------------------------------------------

  parms$chains     <- 5 # number of parallel chains for convergence check
  parms$iters      <- 500
  parms$burnin_imp <- 0#10 # how many imputation iterations should be discarded
  parms$ndt        <- 1#10 # number of imputed datasets to pool esitmaes from (10)
  parms$thin       <- (parms$iters - parms$burnin_imp)/parms$ndt
  parms$pos_dt  <- (parms$burnin_imp+1):parms$iters # candidate datasets (after convergence)
  parms$keep_dt <- parms$pos_dt[seq(1, length(parms$pos_dt), parms$thin)] # keep 1 dataset every thin
  imp_blasso <- impute_BLAS_hans(Z = Xy_mis, 
                                 O = as.data.frame(O),
                                 parms = parms)
  # Store/load results
  saveRDS(object = list(imp_obj = imp_blasso, 
                        parms = parms,
                        cond = cond_convcheck),
          file = paste0(path_loc, "blasso_cnvg", ".rds") )
  
  cnvg_blasso <- readRDS(paste0(path_loc, "blasso_cnvg", ".rds"))
  parms_temp <- cnvg_blasso$parms
  imp_blasso  <- cnvg_blasso$imp_obj
  
  # Convergence Plot
  # Full range
  par(mfrow = c(2,3))
  for (v in 1:length(parms_temp$z_m_id)) {
    imps_4plot <- imp_blasso$imps[[1]][[v]]
    mean_imp <- rowMeans(imps_4plot)
    plot(1:parms_temp$iters, mean_imp, type = "l",
         main = "Blasso Mean Imputations",
         ylim = c(mean(mean_imp)-4*sd(mean_imp), mean(mean_imp)+4*sd(mean_imp)),
         ylab = paste0("z", v), xlab = "Iteration")
    # And add imputations from other chains
    for (i in 2:(parms_temp$chains)) {
      imps_4plot <- imp_blasso$imps[[i]][[v]]
      mean_imp <- rowMeans(imps_4plot)
      lines(1:parms_temp$iters, mean_imp)
    }
  }
  # 0 to 100 iterations
  par(mfrow = c(2,3))
  for (v in 1:length(parms_temp$z_m_id)) {
    imps_4plot <- imp_blasso$imps[[1]][[v]][1:100, ]
    mean_imp <- rowMeans(imps_4plot)
    plot(seq(parms_temp$iters)[1:100], mean_imp, type = "l",
         main = "Blasso Mean Imputations",
         ylim = c(mean(mean_imp)-4*sd(mean_imp), mean(mean_imp)+4*sd(mean_imp)),
         ylab = paste0("z", v), xlab = "Iteration")
    # And add imputations from other chains
    for (i in 2:(parms_temp$chains)) {
      imps_4plot <- imp_blasso$imps[[i]][[v]][1:100, ]
      mean_imp <- rowMeans(imps_4plot)
      lines(seq(parms_temp$iters)[1:100], mean_imp)
    }
  }
  
  # CONCLUSION ?:
  # 20 burn in, 40 iters should be fine in the most extreme case
  # keep 10 dts
  
## Random forest
  parms$chains     <- 5 # number of parallel chains for convergence check
  parms$iters      <- 20
  
  imp_RANF <- impute_RANF(Z = Xy_mis,
                          O = as.data.frame(O),
                          perform = TRUE,
                          parms = parms)
  
  imp_MI_RF_mids <- mice::mice(Xy_mis, 
                               m = parms$iters,
                               maxit = parms$chains,
                               meth = "rf", 
                               ntree = parms$rfntree)
  
  n_iter_plot <- c(parms$iters, 100)[1]
  par(mfrow = c(2,3))
  for (v in 1:length(parms$z_m_id)) {
    imps_4plot <- imp_RANF$imps[[1]][[v]][1:n_iter_plot, ]
    mean_imp <- rowMeans(imps_4plot)
    plot(seq(parms$iters)[1:n_iter_plot], mean_imp, type = "l",
         main = "Blasso Mean Imputations",
         ylim = c(mean(mean_imp)-4*sd(mean_imp), mean(mean_imp)+4*sd(mean_imp)),
         ylab = paste0("z", v), xlab = "Iteration")
    # And add imputations from other chains
    for (i in 2:(parms$chains)) {
      imps_4plot <- imp_RANF$imps[[i]][[v]][1:n_iter_plot, ]
      mean_imp <- rowMeans(imps_4plot)
      lines(seq(parms$iters)[1:n_iter_plot], mean_imp)
    }
  }
  
## impute_PCA
  parms$mice_iters      <- 50
  parms$mice_ndt        <- 5 
  imp_PCA <- impute_PCA(Z = Xy_mis, parms = parms)
  length(imp_PCA$dats)
  plot(imp_PCA$mids)

## imp_MICE_TR
  parms$mice_iters      <- 50
  parms$mice_ndt        <- 5 
  imp_MICE_TR <- impute_MICE_TR(Z = Xy_mis,
                                cond = cond_convcheck,
                                parms = parms)
  length(imp_MICE_TR$dats)
  plot(imp_MICE_TR$mids)
  
  # After 10 iterations everything goes
  