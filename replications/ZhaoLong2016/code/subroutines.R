### Title:    Replication Zhao Long 2016 - Functions
### Author:   Edoardo Costantini
### Created:  2020-02-20
### Modified: 2020-02-20

# Functions ---------------------------------------------------------------


runCell <- function(cond, m = 5, iters = 1) {
  ## Description
  # Generates data for 1 condition and performs imutations accoding to
  # DURR, IURR and BLasso methods.
  ## For internals
  # source("./init.R")
  # source("./functions.R")
  # source("./fun_DURR_impute.R")
  # source("./fun_IURR_impute.R")
  # m = 5
  # iters = 1
  # cond <- conds[1,]

  ## Set up seed ----------------------------------------------------------- ##
  # Setup the PRNG (needs rlecuyer package):
  
  .lec.SetPackageSeed(rep(parms$seed, 6))
  if(!rp %in% .lec.GetStreams())
    .lec.CreateStream(c(1 : parms$nStreams))
  .lec.CurrentStream(rp)
  
  ## Data ------------------------------------------------------------------ ##
  # Gen 1 dataset w/ missing values; 1 fully-obs data for out-of-sample rmse
  
  Xy <- genData(cond, parms)
  Xy_mis <- imposeMiss(Xy, parms)$Xy_miss
  
  miss_descrps <- mean(is.na(Xy_mis[, "z1"]))
  
  ## Imputation ------------------------------------------------------------ ##
  # Impute m times the data w/ missing values w/ different methods
  
  # Impute according to DURR method
  imp_DURR_lasso <- impute_DURR(Xy_mis = Xy_mis, 
                                cond = cond, 
                                chains = m, 
                                iters = iters, 
                                reg_type="lasso")
  
  # MICE default
  imp_IURR_lasso <- impute_DURR(Xy_mis = Xy_mis, 
                                cond = cond, 
                                chains = m, 
                                iters = iters, 
                                reg_type="lasso")
  
  # MICE cart
  imp_BLasso <- impute_DURR(Xy_mis = Xy_mis, 
                                cond = cond, 
                                chains = m, 
                                iters = iters, 
                                reg_type="lasso")
  
  ## Pooling --------------------------------------------------------------- ##
  # Analyse and pool estimates of the m imputations for each method
  
  # DURR
  pool_bs_DURR <- get_pool_est(imp_DURR_lasso, parms)$est
  pool_CI_DURR <- get_pool_est(imp_DURR_lasso, parms)$CI
  
  # IURR
  pool_bs_IURR <- get_pool_est(imp_DURR_lasso, parms)$est # will change to IURR
  pool_CI_IURR <- get_pool_est(imp_DURR_lasso, parms)$CI
  
  # BLasso
  pool_bs_BLas <- get_pool_est(imp_DURR_lasso, parms)$est # will change to BLas
  pool_CI_BLas <- get_pool_est(imp_DURR_lasso, parms)$CI
  
  # collect estiamtes
  pool_est <- cbind(pool_bs_DURR, pool_bs_IURR, pool_bs_BLas)
  pool_conf <- data.frame(DURR = pool_CI_DURR, 
                          IURR = pool_CI_IURR, 
                          BLAs = pool_CI_BLas)
  
  ## Store output ---------------------------------------------------------- ##
  
  output <- list(pool_est = pool_est,
                 pool_conf = pool_conf,
                 miss_descrps = miss_descrps)
  
  return(output)
}
