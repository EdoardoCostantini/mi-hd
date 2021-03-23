### Title:    Exp 5 Rdige Parameter Cross-validation
### Project:  Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2021-03-17

## Make sure we have a clean environment:
rm(list = ls(all = TRUE))

## Initialize the environment:
source("./exp5_init_cv_bridge.R")

## Extract commandline arguments
args <- commandArgs(trailingOnly = TRUE)
rp <- as.numeric(args[1]) # replication rp = 1 to desired
parms$outDir <- args[2]   # overwrite output directory defined in exp5_init.R

## Start Timer
time_1 <- Sys.time()

## Run one replication of the simulation:
  doRep_cluster(rp = rp, conds = conds, parms = parms)
  
## End Timer
time_2 <- Sys.time()
print(time_2 - time_1)
