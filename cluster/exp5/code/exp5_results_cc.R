### Title:    Conergence Check Blasso
### Project:  Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-07-03
### Notes:    The goal is to find out how many iterations should be used
###           in the full study for each method.

rm(list = ls())
source("./init_general.R")

# Input files
resDir <- "../output/exp5_cc/7523352/" # location of results
fileNames <- grep(".rds", list.files(resDir), value = TRUE)
out_new <- lapply(paste0(resDir, fileNames), readRDS)

# Give unique name to all objects
names(out_new) <- fileNames

# Create indexing object 
# Unique repetitions
rep_index <- paste0("rep", 1:10, "[^0-9]")

# Create an index based on the reptition membership
index <- lapply(rep_index, 
                function(x) {
                  grep(x, names(out_new), value = TRUE)
                }
)

# Reaggregate Results 
list_out <- lapply(index, function(x) {
  x <- index[[1]]
  temp_out <- out_new[x]
  names(temp_out) <- gsub("(.*?)cond", "", (names(temp_out))) 
  return(temp_out)
  }
)

# Obtained best Ridge per condition ---------------------------------------
out_cnv <- list_out
out_cnv[[1]]$`lv10pm03fllowridge1e-07.rds`$imp
  
# What to show
exp_dat <- 3 # which data replication (10 possibilities)
iters_range <- 1:250 # which set of iterations
y_range <- c(-1, .5)

# DURR_la (50-100 good for all)
exp_dat <- 2
mean_traceplot(out_cnv, 
               method = "blasso", 
               dat = exp_dat, 
               y_range = y_range, 
               iters = iters_range)

conds_bridge <- bridge_cv(out, 
                          nconds = 4, 
                          nreps = 10, 
                          ridge_values = 10^c(-1, -7))
conds_bridge
