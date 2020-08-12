### Title:    Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-05-19

rm(list=ls())
source("./init_general.R")
source("./exp1_init.R")

## Create a cluster object:
clus <- makeCluster(10)

## Two different ways to source a script on the worker nodes:
clusterEvalQ(cl = clus, expr = source("./init_general.R"))
clusterEvalQ(cl = clus, expr = source("./exp1_init.R"))

## Data directory for storage

# Progress report file ----------------------------------------------------
file.create(paste0(parms$outDir, parms$report_file_name))

cat(paste0("SIMULATION PROGRESS REPORT",
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

# Save output -------------------------------------------------------------

saveRDS(out,
        paste0(parms$outDir,
               parms$results_file_name)
)