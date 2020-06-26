### Title:    Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-05-19

library(parallel) # detectCores(); makeCluster()
rm(list=ls())
source("./init.R")
rp <- 1
i <- 1
## Sequential run
set.seed(12345)

out1 <- list()

for(rp in 1 : 5){
  print(paste0("REPETITION:", rp))
  out1[[rp]] <- doRep(rp, conds = conds, parms = parms)
}

out1[[3]]$cond_500_0.3
out <- out1