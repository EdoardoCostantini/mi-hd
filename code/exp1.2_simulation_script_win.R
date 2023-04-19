# Project:   imputeHD-comp
# Objective: Simulation script for nonlinearity study
# Author:    Edoardo Costantini
# Created:   2023-03-28
# Modified:  2023-04-19
# Notes: 

# Clean environment
rm(list=ls())

# Source scripts
source("./init_general.R")
source("./exp1.2_init.R")

# Create a cluster object:
clus <- makeCluster(30)

# Source the scripts on the worker nodes:
clusterEvalQ(cl = clus, expr = source("./init_general.R"))
clusterEvalQ(cl = clus, expr = source("./exp1.2_init.R"))

# Progress report file ---------------------------------------------------------

# Create an empty file
file.create(paste0(parms$outDir, parms$report_file_name))

# Add headers to the report file
cat(
       paste0(
              "SIMULATION PROGRESS REPORT",
              # Description
              "\n",
              "Description: ", parms$description,
              "\n",
              # Time
              "\n", "------", "\n",
              "Starts at: ", Sys.time(),
              "\n", "------", "\n"
       ),
       file = paste0(parms$outDir, parms$report_file_name),
       sep = "\n",
       append = TRUE
)

# Simulate in parallel ---------------------------------------------------------

# Store time at beginning of the simulation
sim_start <- Sys.time()

# Run the computations in parallel on the 'clus' object:
out <- parLapply(
       cl = clus,
       X = 1:parms$dt_rep,
       fun = doRep,
       conds = conds,
       parms = parms,
       debug = FALSE,
       verbose = FALSE
)

# Kill the cluster:
stopCluster(clus)

# Store time at beginning of the simulation
sim_ends <- Sys.time()

# Attach Info Objects
out$parms <- parms
out$conds <- conds
out$session_info <- devtools::session_info()

# Close report file
cat(
       paste0(
              "\n", "------", "\n",
              "Ends at: ", Sys.time(), "\n",
              "Run time: ",
              round(difftime(sim_ends, sim_start, units = "hours"), 3), " h",
              "\n", "------", "\n"
       ),
       file = paste0(parms$outDir, parms$report_file_name),
       sep = "\n",
       append = TRUE
)

# Save output -------------------------------------------------------------

saveRDS(
       out,
       paste0(
              parms$outDir,
              parms$results_file_name
       )
)