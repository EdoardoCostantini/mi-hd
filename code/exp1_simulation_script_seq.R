### Title:    Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-05-19

library(parallel) # detectCores(); makeCluster()
rm(list=ls())
source("./init_general.R")
source("./exp1_init.R")

rp <- 1
i <- 1
## Sequential run

out1 <- list()

for(rp in 1 : parms$dt_rep){
  print(paste0("REPETITION:", rp))
  out1[[rp]] <- doRep(rp, conds = conds, parms = parms)
}

out$parms <- parms
out$parms$dt_rep <- 3
out <- out1