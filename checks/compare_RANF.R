### Title:    Compare Random Forest Imputation approaches
### Project:  Imputing High Dimensional Data
### Author:   Anonymized for peer review
### Created:  2020-06-26
### Notes:    This shows the equivalence of my fun_RANF_impute.R with an
###           the mice.rf.impute function from mice package.
###           Bias: You can see in the results same bias between the two version
###                 of the imptuation with random forest. 
###           Time: The version I want to use takes less time to impute.
###           Convergence: good, see bottom section

  rm(list=ls())
  source("./init.R")

# Repeated runs -----------------------------------------------------------

  # Modify usual initial script
  parms$dt_rep     <- 100
  
  parms$chains     <- 1 # number of parallel chains == number of imputations
  parms$iters      <- 50
  parms$methods    <- c("RANF", "MICE-RF", "MICE-TR", "GS", "CC")
  parms$ndt        <- 5  # number of imputed datasets to pool esitmaes from (10)
  parms$burnin_imp <- 10 #10 # how many imputation iterations should be discarded
  parms$thin       <- (parms$iters - parms$burnin_imp)/parms$ndt
  parms$pos_dt     <- (parms$burnin_imp+1):parms$iters # candidate datasets (after convergence)
  parms$keep_dt    <- parms$pos_dt[seq(1, length(parms$pos_dt), parms$thin)] # keep 1 dataset every thin
  
  # For mice-like algorithms
  parms$mice_iters      <- 20
  parms$mice_ndt        <- parms$ndt
  
  # Methods
  parms$meth_sel$MI_RF <- TRUE
  
  # Special condition
  p   <- c(50) # number of variables
  rho <- c(0.5) # autoregressive structure
  q   <- c(1) # c(4, 20)  # active set (num of variables are true predictors of y)
  latent <- c(FALSE, TRUE)[1]
  pm <- c(.3)
  
  conds <- expand.grid(p, rho, q, latent, pm)
    colnames(conds) <- c("p", "rho", "q", "latent", "pm")
  
  runCell <- function(cond, parms, rep_status) {
    
    # Generate data
    # cond <- conds[1, ]; set.seed(1234)
    Xy <- genDt_mvn(cond, parms)
    Xy_mis <- imposeMiss(Xy, parms, cond)
    Xy_mis <- cbind(Xy_mis[parms$z_m_id], Xy_mis[-parms$z_m_id])
    
    O <- !is.na(Xy_mis) # matrix index of observed values
    miss_descrps <- colMeans(!O[, 1:parms$zm_n]) 
    
    ## Imputation ------------------------------------------------------------ ##
    
    # Impute using norm and true active set -----------------------------------
    
    imp_MICE_TR <- impute_MICE_TR(Z = Xy_mis,
                                  cond = cond,
                                  parms = parms)
    
    # mice ranf -------------------------------------------------------
    
    imp_MICE_RF <- impute_MICE_RF(Z = Xy_mis,
                                  parms = parms)
    
    # RANF mine ------------------------------------------------------
    
    imp_RANF <- impute_RANF(Z = Xy_mis,
                            O = as.data.frame(O),
                            perform = TRUE,
                            parms = parms)
    
    ## Analyse --------------------------------------------------------------- ##
    # For each imp method, analyse all datasets based on model defined in init.R
    fits_md <- lapply(list(imp_RANF$dats,
                           imp_MICE_RF$dats,
                           imp_MICE_TR$dats), 
                      fit_sat_model)
    
    # Gold Standard
    fit_gs <- sem(parms$lav_model, 
                  data = Xy, 
                  likelihood = "wishart")
    
    # CC (complete case analysis)
    fit_cc <- sem(parms$lav_model, 
                  data = Xy_mis[rowSums(!O) == 0, ], # select complete casese
                  likelihood = "wishart")
    
    ## Pooling --------------------------------------------------------------- ##
    # For each imp method, pool estimates across the m datasets
    pool_EST <- sapply(fits_md[lapply(fits_md, length) != 0], 
                       get_pool_EST)
    pool_CI  <- sapply(fits_md[lapply(fits_md, length) != 0], 
                       get_pool_CI)
    
    # append single imputations, and GS and CC results
    all_EST <- cbind(pool_EST, 
                     get_EST(list(fit_gs, fit_cc)))
    
    all_CI <- cbind(pool_CI, 
                    get_CI(list(fit_gs, fit_cc)))
    
    # give meaningful names
    colnames(all_EST) <- parms$methods
    colnames(all_CI) <- parms$methods
    
    # compute bias and check coverage
    # cond_bias <- sapply(all_EST, bias_est, x_true = parms$b)
    # cond_CIco <- sapply(pool_CI, check_cover)
    
    # aggregate times
    imp_time <- list(ranf_imp    = imp_RANF$time,
                     MICE_RF     = imp_MICE_RF$time,
                     MICE_TR     = imp_MICE_TR$time)
    imp_time <- do.call(cbind, imp_time)[1,]
    
    ## Store output ---------------------------------------------------------- ##
    output <- list(cond = cond,
                   dat_full = Xy,
                   dat_miss = Xy_mis,
                   all_EST = all_EST,
                   all_CI = all_CI,
                   miss_descrps = miss_descrps,
                   run_time_min = imp_time)
    return(output) 
  }
  
  doRep <- function(rp, conds, parms) { 
    ## Setup the PRNG: (for setting seed with guaranteed independece)
    # form the rlecuyer package
    # rp is the replication index, in the mcapply function it the argument
    # that takes the elements of the vector specified for X
    # This creates names for the streams so that you can recover them
    .lec.SetPackageSeed(rep(parms$seed, 6))
    .lec.CreateStream(c(1 : parms$nStreams)) # create 1000 streams
    .lec.CurrentStream(rp) # use the rp sequence out of the 1000
    
    ## Do the calculations for each set of crossed conditions:
    rp_out <- vector("list", nrow(conds)) # store output of repetition
    names(rp_out) <- paste0("cond_", paste0(conds[, 1], "_",  conds[, 3]) )
    
    for(i in 1 : nrow(conds)) {
      exit_while <- "no"
      
      while(exit_while == "no"){
        
        rp_out[[i]] <- try(runCell(cond = conds[i, ], 
                                 parms = parms,
                                 rep_status = rp),
                           silent = TRUE)
        
        if (class(rp_out[[i]]) != "try-error") {
          exit_while <- "yes"
        }
    }
    }
    
    ## Return Function output  
    return(rp_out)
  }
  
  # Run simulation
  out <- mclapply(X        = 1 : parms$dt_rep,
                  FUN      = doRep,
                  conds    = conds,
                  parms    = parms,
                  mc.cores = ( 10 ) )
  out$parms <- parms
  
  # Save/load results
  saveRDS(out, "../checks/output/compare_RANF.rds")
  out <- readRDS("../checks/output/compare_RANF.rds")
  
  # Results
  select_cond <- names(out[[2]])[1] # available conditions
  
  # Step 1. Obtain "true" comparison values
  
  full_dat_est <- NULL
  
  for (i in 1:out$parms$dt_rep) {
    dat <- out[[i]][[select_cond]]$dat_full
    fit <- sem(parms$lav_model, 
               data = dat,
               likelihood = "wishart")
    full_dat_est <- rbind(full_dat_est, 
                          parameterEstimates(fit)$est)
  }
  
  psd_tr_vec <- colMeans(full_dat_est) # pseudo true values
  
  # Step 2. Compute averages of statistics
  
  # Row index for type of paramter
  avg_indx <- 1:out$parms$zm_n
  var_indx <- (out$parms$zm_n+1):(out$parms$zm_n*2)
  cov_indx <- (tail(var_indx, 1)+1):nrow(out[[i]][[select_cond]]$all_EST)
  
  # Store objects
  sum_stats <- matrix(0, 
                      nrow = nrow(out[[i]][[select_cond]]$all_EST), 
                      ncol = length(out$parms$methods))
  
  # Compute averages of the statistics
  for (i in 1:out$parms$dt_rep) {
    sum_stats <- sum_stats + out[[i]][[select_cond]]$all_EST
  }
  
  avg_stats <- sum_stats / out$parms$dt_rep
  
  # Step 3. Obtain Bias
  # Raw bias
  bias <- avg_stats - psd_tr_vec
  round(bias, 3)
  
  # Bias as percentage of true value
  round(bias/psd_tr_vec*100, 1)
  
  # Time
  res_time <- NULL
  for (i in 1:out$parms$dt_rep) {
    if(class(out[[i]]) != "try-error"){
      res_time <- rbind(res_time, out[[i]][[select_cond]]$run_time_min)
    }
  }
  round(colMeans(res_time), 3) # takes half the time or less
  
# Convergence checks ------------------------------------------------------

  # Gen data
  Xy <- genDt_mvn(cond[1,], parms)
  Xy_mis <- imposeMiss(Xy, parms, cond)
  Xy_mis <- cbind(Xy_mis[parms$z_m_id], Xy_mis[-parms$z_m_id])
  
  O <- !is.na(Xy_mis) # matrix index of observed values
  miss_descrps <- colMeans(!O[, 1:parms$zm_n]) 
  
  
  # Run imputations
  imp_MI_RF_mids <- mice::mice(Xy_mis, 
                               m = parms$mice_ndt,
                               maxit = parms$mice_iters,
                               meth = "rf", 
                               ntree = parms$rfntree)
  imp_RANF <- impute_RANF(Z = Xy_mis,
                          O = as.data.frame(O),
                          perform = TRUE,
                          parms = parms)
  
  # Convergence
  # Regular mice
  plot(imp_MI_RF_mids)
  
  # Your Ranf
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
  