### Title:    Conergence Check Blasso
### Project:  Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-07-03
### Notes:    The goal is to find out how many iterations should be used
###           in the full study for each method.

rm(list = ls())
source("./init_general.R")

# Input files
resDir <- "../output/exp5_cv_bridge_7523347/7523347/" # location of results
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
conds_bridge <- bridge_cv(out, 
                          nconds = 4, 
                          nreps = 10, 
                          ridge_values = 10^c(-1, -7))
conds_bridge
