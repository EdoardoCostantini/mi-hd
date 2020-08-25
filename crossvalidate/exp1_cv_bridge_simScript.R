### Title:    Simulation script for ridge paramter cross-validation
### Project:  Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-08-25
### Notes:    The goal is to find out how many iterations should be used
###           in the full study for each method.

rm(list = ls())
source("./init_general.R")
source("./exp1_init.R")
source("../crossvalidate/exp1_cv_bridge_init.R")

sim_start <- Sys.time()

# Create a cluster object:
clus <- makeCluster(parms$dt_rep)

# Prep Environments
clusterEvalQ(cl = clus, expr = source("./init_general.R"))
clusterEvalQ(cl = clus, expr = source("./exp1_init.R"))
clusterEvalQ(cl = clus, expr = source("../crossvalidate/exp1_cv_bridge_init.R"))

## Run the computations in parallel on the 'clus' object:
out <- parLapply(cl = clus, 
                 X = 1 : parms$dt_rep,
                 fun = doRep, 
                 conds = conds, 
                 parms = parms)

## Kill the cluster:
stopCluster(clus)

sim_ends <- Sys.time()

## Append Parms and conds to out object
out$parms <- parms
out$conds <- conds

## Save results
saveRDS(out,
        paste0(parms$outDir,
               parms$results_file_name)
)