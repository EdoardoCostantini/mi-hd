### Title:    Simulation script for ridge paramter cross-validation
### Project:  Imputing High Dimensional Data
### Author:   Anonymized for peer review
### Created:  2020-08-24
### Notes:    The goal is to find out how many iterations should be used
###           in the full study for each method.

rm(list = ls())
source("./init_general.R")
source("./exp4_init.R")
source("../crossvalidate/exp4_cv_bridge_init.R")

# Create a cluster object:
clus <- makeCluster(parms$dt_rep)

# Prep Environments
clusterEvalQ(cl = clus, expr = source("./init_general.R"))
clusterEvalQ(cl = clus, expr = source("./exp4_init.R"))
clusterEvalQ(cl = clus, expr = source("../crossvalidate/exp4_cv_bridge_init.R"))

# Progress report file ----------------------------------------------------
file.create(paste0(parms$outDir, parms$report_file_name))

cat(paste0("BRIDGE CROSSVALIDATION PROGRESS REPORT",
           "\n",
           "Description: ", parms$description, "\n",
           "\n", "------", "\n",
           "Starts at: ", Sys.time(),
           "\n", "------", "\n" ),
    file = paste0(parms$outDir, parms$report_file_name),
    sep = "\n",
    append = TRUE)

## Start Time
sim_start <- Sys.time()

## Run the computations in parallel on the 'clus' object:
out <- parLapply(cl = clus, 
                 X = 1 : parms$dt_rep,
                 fun = doRep, 
                 conds = conds, 
                 parms = parms)

## Kill the cluster:
stopCluster(clus)

sim_ends <- Sys.time()

cat(paste0("\n", "------", "\n",
           "Ends at: ", Sys.time(), "\n",
           "Run time: ",
           round(difftime(sim_ends, sim_start, units = "hours"), 3), " h",
           "\n", "------", "\n"),
    file = paste0(parms$outDir, parms$report_file_name),
    sep = "\n",
    append = TRUE)

## Append Parms and conds to out object
out$parms <- parms
out$conds <- conds

## Save results
saveRDS(out,
        paste0(parms$outDir,
               parms$results_file_name)
)