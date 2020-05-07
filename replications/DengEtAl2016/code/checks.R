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
  for (i in 1:parms$dt_rep) {
    Xy <- genData(cond, parms)
    Xy_mis <- imposeMiss(Xy, parms)$Xy_miss
    obj[i] <- mean(rowSums(is.na(Xy_mis)) != 0)
  }
  round(mean(obj), 1) # approximately .4
  # Should be .4 missingness on this variable on average