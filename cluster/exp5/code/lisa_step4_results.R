### Title:    Conergence Check Blasso
### Project:  Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-07-03
### Notes:    The goal is to find out how many iterations should be used
###           in the full study for each method.

  rm(list = ls())
  source("./init_general.R")

# Load results ------------------------------------------------------------

  # Output directoru
  outDir <- "../output/"
  
  # Job ID
  idJob <- c("7652179", "7653072")[2]
  
  # Define Output Run folder
  resDir <- paste0(outDir, "exp5_mainSim_", idJob, "/") # location of results
  resDir <- paste0(outDir, idJob, "/") # location of results
  
  # Create Empty Directory to contain every output
  subDir_tar <- grep(".tar.gz", list.files(resDir), value = TRUE)
  subDir     <- paste0(resDir, "rds_objects/")
  system(command = paste0("mkdir -p ", subDir))
  
  # Untar everything and move to main output folder
  untar_commands <- paste0("tar ", 
                           "-C ", 
                           subDir, # where to place
                           " -xvf ", # untarring
                           resDir, subDir_tar)
  lapply(untar_commands, system)
  
  # Obtain unique names of all result files
  fileNames <- grep(".rds", list.files(subDir), value = TRUE)
  
  # Read all of them in
  out <- lapply(paste0(subDir, fileNames), readRDS)
  
  # Read the session info object
  sInfo <- readRDS(paste0(resDir, "sInfo.rds"))
  
# Restructure Results -----------------------------------------------------

  # Give unique name to all objects
  names(out) <- fileNames
  
  # Define the Unique repetitions
  reps <- sInfo$parms$dt_rep # you need to input this manually right now
  rep_index <- paste0("rep", 1:reps, "[^0-9]")
  
  # Create an index based on the reptition membership
  index <- lapply(rep_index, 
                  function(x) {
                    grep(x, names(out), value = TRUE)
                  }
  )
  
  # Reaggregate Results 
  out_list <- lapply(index, function(x) {
    temp_out <- out[x]
    names(temp_out) <- gsub("(.*?)cond", "", (names(temp_out))) 
    return(temp_out)
    }
  )
  out <- out_list
  
  # Append the parms object for this run
  out$parms <- sInfo$parms
  out$conds <- sInfo$conds
  out$session_info <- sInfo$session_info

# Time --------------------------------------------------------------------

  out_time <- sapply(1:length(names(out[[1]])), res_sem_time, out = out)
  colnames(out_time) <- names(out[[1]])
  t(out_time)

# Extract Estimates ------------------------------------------------------
# Extract results per conditions for each analysis type

## SEM estiamtes raw data (saturated model) ##
  semR_res <- lapply(seq_along(1:length(out[[1]])),
                     function(x) res_sum(out, 
                                         model = "semR", 
                                         condition = x,
                                         bias_sd = TRUE))

## CFA model results
  CFA_res <- lapply(1:length(out[[1]]),
                    function(x) res_sum(out, 
                                        model = "CFA", 
                                        condition = x))

  # Check
  lapply(1:length(out[[1]]),
         function(x) CFA_res[[x]]$bias_per[1:10, ])

## SEM estaimted Scored data
  semS_res <- lapply(1:length(out[[1]]),
                     function(x) res_sum(out, 
                                         model = "semS", 
                                         condition = x,
                                         bias_sd = TRUE))

  # Show results for a given condition
  lapply(1:length(out[[1]]),
         function(x) semS_res[[x]]$bias_raw[1:2,])

# Save Results ------------------------------------------------------------

  output <- lapply(list(semR = semR_res,
                        CFA  = CFA_res,
                        semS = semS_res), 
                   function(x){
                     names(x) <- paste0("cond", seq_along(out[[1]]))
                     return(x)
                   }
  )
  output <- c(output, parms = list(out$parms), conds = list(out$conds))
  
  
  