### Title:    Imputing High Dimensional Data: Addendum
### Author:   Edoardo Costantini
### Created:  2021-01-29

## Make sure we have a clean environment:
rm(list = ls(all = TRUE))

## Initialize the environment:
source("./init_general.R")
source("./exp5_init.R")

## Extract commandline arguments
args <- commandArgs(trailingOnly = TRUE)
rp <- as.numeric(args[1]) # replication rp = 1 to desired
parms$outDir <- args[2]   # overwrite output directory defined in exp5_init.R

## Run one replication of the simulation:
doRep_cluster(rp = rp, conds = conds, parms = parms)

## Time Estiamte
start <- Sys.time()
out <- doRep(rp = 1, conds = conds, parms = parms,
             verbose = TRUE)
end <- Sys.time()
end-start
sink() # end sing