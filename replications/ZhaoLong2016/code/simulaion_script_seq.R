### Title:    Replication Zhao Long 2016 - Simulation script
### Author:   Anonymized for peer review
### Created:  2020-03-03
### Modified: 2020-03-03

library(parallel) # detectCores(); makeCluster()

source("./init.R")

## Sequential run
out1 <- list()
for(rp in 1 : parms$dt_rep){
  out1[[rp]] <- doRep(rp, conds = conds, parms = parms)
}
