### Title:    Replication Zhao Long 2016 - Simulation script
### Author:   Anonymized for peer review
### Created:  2020-03-03
### Modified: 2020-03-03

library(parallel) # detectCores(); makeCluster()
rm(list=ls())
source("./init.R")

## Sequential run
set.seed(1234)
out1 <- list()
for(rp in 1 : parms$dt_rep){
  print(paste0("REPETITION:", rp))
  out1[[rp]] <- try(doRep(rp, conds = conds, parms = parms), silent = TRUE)
}
