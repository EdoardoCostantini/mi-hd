### Title:    Replication Howard et al 2015
### Author:   Anonymized for peer review
### Created:  2020-05-18

# Functions ---------------------------------------------------------------

## Run one replication of the simulation:
doRep <- function(rp, conds, parms) { 
  ## Setup the PRNG: (for setting seed with guaranteed independece)
  # form the rlecuyer package
  # rp is the replication index, in the mcapply function it the argument
  # that takes the elements of the vector specified for X
  # This creates names for the streams so that you can recover them
  .lec.SetPackageSeed(rep(parms$seed, 6))
  .lec.CreateStream(c(1 : parms$nStreams)) # create 1000 streams
  .lec.CurrentStream(rp) # use the rp sequence out of the 1000
  
  ## Progress report - Start
  cat(paste0(Sys.time(), " - Starts Repetition: ", rp, 
             "\n"),
      file = paste0(parms$outDir, parms$report_file_name),
      append = TRUE)
  
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

runCell <- function(cond, parms, rep_status) {
  ## Description
  # Generates 1 data for 1 condition and estiamtes
  # descriptive statistics (mean, variance, covarainces)
  # using FIML
  ## For internals
  # cond <- conds[1, ]
  
  ## Data ------------------------------------------------------------------ ##
  # Gen 1 population dataset dataset and 1 sampled version of it
  Xy <- genData(cond, parms)
  
  # Impose missingness
  if(cond["misMch"] == 1L){
    Xy_mis <- impose_L_Miss(Xy$smp, cond, parms)
  } else {
    Xy_mis <- impose_NL_Miss(Xy$smp, cond, parms)
  }
  
  Xy_mis <- impose_MCAR(Xy_mis, cond, parms)
  
  O <- !is.na(Xy_mis) # matrix index of observed values
  miss_descrps <- colMeans(!O) # check correct miss %
  
  ## Gold Standards ##
  fit_GSp <- sem(mod_IAV, Xy$pop)
  pe_GSp <- parameterEstimates(fit_GSp)[,-c(5, 6, 7)]
  
  #
  cat(paste0(Sys.time(), " | Reptetition ", rep_status, ": GSp done",
             "\n"),
      file = paste0(parms$outDir, parms$report_file_name),
      append = TRUE)
  
  fit_GSs <- sem(mod_IAV, Xy$smp)
  pe_GSs <- parameterEstimates(fit_GSs)[,-c(5, 6, 7)]
  
  #
  cat(paste0(Sys.time(), " | Reptetition ", rep_status, ": GSs done",
             "\n"),
      file = paste0(parms$outDir, parms$report_file_name),
      append = TRUE)
  
  ## FIML Analysis ##
  # No Auxiliary variables
  fit_NAV <- sem(mod_NAV, Xy_mis, missing = 'fiml')
  pe_NAV <- parameterEstimates(fit_NAV)[,-c(5, 6, 7)]
  
  cat(paste0(Sys.time(), " | Reptetition ", rep_status, ": NAV done",
             "\n"),
      file = paste0(parms$outDir, parms$report_file_name),
      append = TRUE)
  
  # Inclusive Auxiliary vairabels
  fit_IAV <- sem(mod_IAV, Xy_mis, missing = 'fiml')
  pe_IAV <- parameterEstimates(fit_IAV)[,-c(5, 6, 7)]
  
  cat(paste0(Sys.time(), " | Reptetition ", rep_status, ": IAV done",
             "\n"),
      file = paste0(parms$outDir, parms$report_file_name),
      append = TRUE)
  
  # PCA approach
  # Fill in values single random prediction
  imp <- mice(Xy_mis, 
              method = "norm.predict", 
              m = 1, maxit = 1,
              printFlag = FALSE)
  imp_dats <- mice::complete(imp, "all")$'1'
  
  # Extract PCA
  pr.out <- prcomp(imp_dats[, -c(1,2)], scale. = TRUE)
  pr.out$sdev^2/sum(pr.out$sdev^2) # check explained variance
  PC1 <- pr.out$x[,1]
  dt_pcax <- cbind(Xy_mis[, 1:2], PC1)
  
  # Fit FIML
  fit_PCA <- sem(mod_PCA, dt_pcax, missing = 'fiml')
  pe_PCA <- parameterEstimates(fit_PCA)[,-c(5, 6, 7)]
  
  cat(paste0(Sys.time(), " | Reptetition ", rep_status, ": PCA done",
             "\n"),
      file = paste0(parms$outDir, parms$report_file_name),
      append = TRUE)
  
  ## Results ##
  # Gather results
  pe <- list(pe_NAV, pe_IAV, pe_PCA, pe_GSp, pe_GSs)
  rep_bias <- sapply(pe, bias, parms = parms)
  rep_cov <- sapply(pe, cove, parms = parms)
  
  colnames(rep_bias) <- parms$methods
  colnames(rep_cov) <- parms$methods
  
  ## Store output ---------------------------------------------------------- ##
  output <- list(rep_bias = rep_bias,
                 rep_cov = rep_cov,
                 parms = parms)
  return(output)
}
