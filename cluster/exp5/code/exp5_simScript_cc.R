### Title:    Exp 5 Convergence Checks (script for lisa)
### Project:  Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2021-03-17


# Lisa Version ------------------------------------------------------------
# Comment out if you want to use the windows version
# 
# ## Make sure we have a clean environment:
# rm(list = ls(all = TRUE))
# 
# ## Initialize the environment:
# source("./exp5_init_cc.R")
# 
# ## Extract commandline arguments
# args <- commandArgs(trailingOnly = TRUE)
# rp <- as.numeric(args[1]) # replication rp = 1 to desired
# parms$outDir <- args[2]   # overwrite output directory defined in exp5_init.R
# 
# ## Start Timer
# time_1 <- Sys.time()
# 
# ## Run one replication of the simulation:
# doRep_cluster(rp = rp, conds = conds, parms = parms)
# 
# ## End Timer
# time_2 <- Sys.time()
# print(time_2 - time_1)

# Windwos -----------------------------------------------------------------
# Comment out if you want to use the Lisa version

rm(list = ls())
source("./init_general.R")
source("./exp5_init.R")
source("./exp5_init_cc.R")

sim_start <- Sys.time()

# Create a cluster object:
clus <- makeCluster(parms$dt_rep)

# Prep Environments
clusterEvalQ(cl = clus, expr = source("./init_general.R"))
clusterEvalQ(cl = clus, expr = source("./exp5_init.R"))
clusterEvalQ(cl = clus, expr = source("./exp5_init_cc.R"))

## Run the computations in parallel on the 'clus' object:
out <- parLapply(cl = clus, 
                 X = 1 : parms$dt_rep,
                 fun = doRep, 
                 conds = conds, 
                 parms = parms)

## Kill the cluster:
stopCluster(clus)

sim_ends <- Sys.time()

out$parms <- parms
out$conds <- conds

## Save results
saveRDS(out,
        paste0(parms$outDir,
               parms$results_file_name)
)
