### Title:    Replication Zhao Long 2016 - Check certain key elements
### Author:   Edoardo Costantini
### Created:  2020-02-20
### Modified: 2020-02-20
### Notes:    Check whether some setups of the simulation are achieved

source("./init.R")
source("./functions.R")

# Missingness -------------------------------------------------------------

  cond <- conds[1,]
  obj <- NULL
  for (i in 1:1e3) {
    Xy <- genData(cond, parms)
    obj[i] <- mean(imposeMiss(Xy, parms)$nR)
  }
  mean(obj)
  # Should be .3 missingness on this variable on average