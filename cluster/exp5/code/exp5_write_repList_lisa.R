### Title:   Create replication list
### Author:  Edoardo Costantini
### Created: 2021-03-16
### Notes:   Needed for Lisa cluster run

# Stopos Pool for Bridge Cross-validation ---------------------------------
  source("./exp5_init_cv_bridge.R")
  
  rep_start <- 1
  rep_end   <- parms$dt_rep
  
  # Save in file for Stopos
  outDir <- "../input/"
  fileName <- paste0("exp", parms$exp, 
                     "_SPL",
                     "_cv_bridge")
  write(as.character(rep_start:rep_end),
        file = paste0(outDir, fileName))
  

# Stopos Pool for Covergence Checks ---------------------------------------
  source("./exp5_init_cc.R")
  
  rep_start <- 1
  rep_end   <- parms$dt_rep
  
  # Save in file for Stopos
  outDir <- "../input/"
  fileName <- paste0("exp", parms$exp, 
                     "_SPL",
                     "_cc")
  write(as.character(rep_start:rep_end),
        file = paste0(outDir, fileName))
  
# Stopos Pool for main simulation -----------------------------------------
  source("./exp5_init.R")
  
  # Calculate
  ncores    <- 15 # I want to use 16 cores in each node
  narray    <- 2  # I want to specify a sbatch array of 2 tasks (sbatch -a 1-2 job_script_array.sh)
  goal_reps <- ncores*narray # should match your total goal of clusters
  rep_start <- 1
  rep_end   <- ncores*narray
  
  # Save in file for Stopos
  outDir <- "../input/"
  fileName <- paste0("exp", parms$exp, 
                     "_SPL", # for Stopos Pool Lines
                     "_simOut")
  write(as.character(rep_start:rep_end),
        file = paste0(outDir, fileName))
  