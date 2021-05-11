### Title:    Imputing High Dimensional Data: Addendum
### Project:  Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2021-03-16
### Notes:    "ms" in title stands for "Main Simulation"

## Make sure we have a clean environment:
rm(list = ls(all = TRUE))

## Initialize the environment:
source("./init_general.R")
source("./exp5_init.R")

## Extract commandline arguments
args <- commandArgs(trailingOnly = TRUE)
rp   <- as.numeric(args[1]) # replication rp = 1 to desired
parms$outDir <- args[2]   # overwrite output directory defined in exp5_init.R

## Run one replication of the simulation:
doRep_cluster(rp = rp, conds = conds, parms = parms)