### Title:   Create replication list
### Author:  Edoardo Costantini
### Created: 2021-03-16
### Notes:   Needed for Lisa cluster run

# Prepare Environment
rm(list = ls())  
source("./exp5_init.R")

# Calculate
ncores    <- 15 # I want to use 15 cores in each node
narray    <- 2  # I want to specify a sbatch array of 2 tasks (sbatch -a 1-2 job_script_array.sh)
goal_reps <- ncores*narray # should match your total goal of repetitions
rep_start <- 1
rep_end   <- goal_reps

# Save in file for Stopos
outDir <- "../input/"
fileName <- paste0("stopos_lines")
write(as.character(rep_start:rep_end),
      file = paste0(outDir, fileName))

# Copmute Estimated CPU time (not printed, just for yourself)
n_nodes <- goal_reps/ncores # this is also the number of nodes
n_cores <- ncores
job_time <- 35 * 1.5 # 50ish hours
n_nodes * n_cores * job_time