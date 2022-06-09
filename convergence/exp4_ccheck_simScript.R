### Title:    Convergence Checks: Initialization Script for Exp 2 (windows)
### Project:  Imputing High Dimensional Data
### Author:   Anonymized for peer review
### Created:  2020-08-11
### Notes:    Runs the simulation for convergence checks (few repetitions, 
###           more chains)
### IMPORTANT: Working directory: ~/imputeHD-comp/code

rm(list = ls())
source("./init_general.R")
source("../convergence/exp4_ccheck_init.R")

## Create a cluster object:
clus <- makeCluster(parms$dt_rep)

## Two different ways to source a script on the worker nodes:
clusterEvalQ(cl = clus, expr = source("./init_general.R"))
clusterEvalQ(cl = clus, expr = source("../convergence/exp4_ccheck_init.R"))

## Data directory for storage

# Progress report file ----------------------------------------------------
file.create(paste0(parms$outDir, parms$report_file_name))

cat(paste0("CONVERGENCE CHECK PROGRESS REPORT",
           "\n",
           "Description: ", parms$description, "\n",
           "\n", "------", "\n",
           "Starts at: ", Sys.time(),
           "\n", "------", "\n" ),
    file = paste0(parms$outDir, parms$report_file_name),
    sep = "\n",
    append = TRUE)

# mcApply parallel --------------------------------------------------------

sim_start <- Sys.time()

## Run the computations in parallel on the 'clus' object:
out <- parLapply(cl = clus, 
                 X = 1 : parms$dt_rep,
                 fun = doRep, 
                 conds = conds, 
                 debug = FALSE,
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

# Attach parm object
out$parms <- parms
out$conds <- conds
out_cnv <- out

# Save output -------------------------------------------------------------

saveRDS(out,
        paste0(parms$outDir,
               parms$results_file_name)
)
