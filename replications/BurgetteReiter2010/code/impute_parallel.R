### Title:    Run imputations (parallel)
### Author:   Edoardo Costantini
### Created:  2020-02-13
### Modified: 2020-02-15
### Note:     This version of the parallel imputation uses parLapply
###           Howver, I was not able to implement the seed. I get the
###           same results out of each repetition. With mcapply no problem

library(parallel) # detectCores(); makeCluster()

source("./functions.R") # applyLib cluster set uo
source("./init.R")      # the parLapply function takes argument from the global
                        # environment

## Data directory for storage

outDir <- "../output/"

# Initialize Environments for parallel ------------------------------------

## Initialize local environments
  # Define environments
  cores <- detectCores() - 1
  cl <- makeCluster(cores)
  
  # Define a vector of necessary packages (to broadcast to local environments)
  packages <- c("mvtnorm", "mice", "miceImpHDv")
  
  # Setup the environment on worker nodes:
  clusterCall(cl = cl, fun = applyLib, pkgList = packages) # broadcast packages
  clusterCall(cl = cl, fun = source, file = "init.R")
  clusterCall(cl = cl, fun = source, file = "functions.R")
  clusterCall(cl = cl, fun = source, file = "simMissingness.R")
  clusterCall(cl = cl, fun = source, file = "subroutines.R")
  
# Run simulation in parallel ----------------------------------------------

## Run simulation
  time1 <- Sys.time()
  out <- parLapply(cl = cl, 
                   X = as.list(1:parms$dt_rep), 
                   fun = function(x) singleRun(chains = parms$chains, 
                                               iters = parms$iters,
                                               parms = parms))
  runTime <- Sys.time() - time1

  stopCluster(cl) # stop the threads

# Save output -------------------------------------------------------------
  
saveRDS(out,
        paste0(outDir,
               "pooled_BR2010-",
               format(Sys.time(), "%Y%m%d_%H%M"),
               ".rds")
)