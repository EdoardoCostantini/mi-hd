### Title:    comparing IURR implementations within mice framework and outside
### Author:   Anonymized for peer review
### Created:  2020-06-00
### Notes:    This shows the equivalence of my fun_IURR_impute.R with an
###           implementation of IURR within the mice R package framework

rm(list=ls())
source("./init.R")
source("../checks/mice.impute.IURR.R") #DURR incorporated in mice

# Modify usual initial script
parms$dt_rep     <- 100
parms$methods    <- c("imp.norm", 
                      "imp.norm.tr", 
                      "IURR", # DengEtAl 2016 strict implementation
                      "imp.IURR", # w/in mice method
                      "GS", "CC") # always there
parms$burnin_imp <- 10 # how many imputation iterations should be discarded
parms$ndt        <- 10 # number of imputed datasets to pool esitmaes from (10)
parms$iters      <- 20
parms$thin       <- (parms$iters - parms$burnin_imp)/parms$ndt
parms$pos_dt     <- (parms$burnin_imp+1):parms$iters # candidate datasets (after convergence)
parms$keep_dt    <- parms$pos_dt[seq(1, length(parms$pos_dt), parms$thin)] # keep 1 dataset every thin

runCell <- function(cond, parms, rep_status) {

# Generate data
  # cond <- conds[1, ]; set.seed(1234)
  Xy <- genData(cond, parms)
  Xy_mis <- imposeMiss(Xy, parms)$Xy_miss
  Xy <- Xy[, c(1:60, 201)]
  Xy_mis <- Xy_mis[, c(1:60, 201)]
  O <- !is.na(Xy_mis) # matrix index of observed values  

## Imputation ------------------------------------------------------------ ##

  # Modify init script for mice like imputations
  parms$chains     <- 10 # number of parallel chains == number of imputations
    
# Impute using norm and true active set -----------------------------------

  imp_MICE_TR <- impute_MICE_TR(Z = Xy_mis,
                                AS_size = cond[3],
                                parms = parms)
  
# Impute using norm -------------------------------------------------------

  imp_MITR_mids <- mice::mice(Xy_mis, 
                              m = parms$ndt,
                              maxit = parms$iters,
                              ridge = 1e-5,
                              method = "norm")
  
# IURR w/in mice pack -----------------------------------------------------
# custom DURR method inside mice package framework
  
  imp_IURR_mids <- mice::mice(Xy_mis,
                              m = parms$ndt,
                              maxit = parms$iters,
                              ridge = 1e-5,
                              method = "IURR")

# DURR (function v1) ------------------------------------------------------
  parms$chains     <- 1 # number of parallel chains == number of imputations
  
  imp_IURR <- impute_IURR(Z = Xy_mis,
                          O = as.data.frame(O),
                          reg_type = "lasso",
                          perform = TRUE,
                          parms = parms)
  
# Compare results ---------------------------------------------------------
  
  fits_md <- lapply(list(mice::complete(imp_MITR_mids, "all"), # imp.norm
                         imp_MICE_TR$dats, # imp.norm.tr
                         imp_IURR$dats, # DURR v1
                         mice::complete(imp_IURR_mids, "all") # imp.DURR 
                         ), 
                    fit_models, mod = parms$formula)
  
  # Gold Standard
  lm_fit_gs <- lm(parms$formula, data = Xy)
  bs_gs <- coef(lm_fit_gs)
  CI_gs <- confint(lm_fit_gs)
  
  # CC (complete case analysis)
  lm_fit_cc <- lm(parms$formula, data = Xy_mis, na.action = na.omit)
  bs_cc <- coef(lm_fit_cc)
  CI_cc <- confint(lm_fit_cc)
  
  ## Pooling --------------------------------------------------------------- ##
  # For each imp method, pool estimates across the m datasets
  pool_EST <- lapply(fits_md[lapply(fits_md, length) != 0], 
                     get_pool_EST)
  pool_CI  <- lapply(fits_md[lapply(fits_md, length) != 0], 
                     get_pool_CI)
  
  
  # append GS and CC results
  pool_EST <- c(pool_EST, list(bs_gs, bs_cc))
  pool_CI <- c(pool_CI, list(CI_gs, CI_cc))
  
  # give meaningful names
  names(pool_EST) <- parms$methods
  names(pool_CI) <- parms$methods
  
  cond_bias <- sapply(pool_EST, bias_est, x_true = parms$b)
  cond_CIco <- sapply(pool_CI, check_cover)
  
  # Return output
  output <- list(pool_est = pool_EST,
                 pool_conf = pool_CI,
                 cond_bias = cond_bias,
                 cond_CIco = cond_CIco,
                 cond = cond,
                 parms = parms)
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
      rp_out[[i]] <- runCell(cond = conds[i, ], 
                                 parms = parms,
                                 rep_status = rp)
    }
  ## Return Function output  
  return(rp_out)
}

# Run simulation
out <- mclapply(X        = 1 : parms$dt_rep,
                FUN      = doRep,
                conds    = conds,
                parms    = parms,
                mc.cores = ( 15 ) )

extract_results(cond_name = names(out[[1]])[1], 
                output = out, 
                dt_rep = out[[1]]$cond_200_4$parms$dt_rep)

# Save results
saveRDS(out,
        paste0(parms$outDir,
               "compare_IURR_100reps_",
               Sys.Date(),
               ".rds")
)

# Load results
out_IURR <- readRDS("../output/compare_IURR_100reps_2020-06-09.rds")

extract_results(cond_name = names(out_IURR[[1]])[1], 
                output = out_IURR, 
                dt_rep = out_IURR[[1]]$cond_200_4$parms$dt_rep)
