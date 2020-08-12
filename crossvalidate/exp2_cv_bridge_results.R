### Title:    Conergence Check Blasso
### Project:  Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-07-03
### Notes:    The goal is to find out how many iterations should be used
###           in the full study for each method.

rm(list = ls())
source("./init_general.R")
source("./exp2_init.R")
source("../crossvalidate/exp2_cv_bridge_init.R")

# Extract average fmi for a given model type ------------------------------
out <- readRDS("../output/exp2_cv_bridge_20200812_1132.rds")

store_0 <- NULL
for (i in 1:nrow(out$conds)) {
  store_1 <- NULL
  for (dt in 1:out$parms$dt_rep) {
    store_1 <- cbind(store_1, c(out[[dt]][[i]]$fmi$semR,
                                #out[[dt]][[i]]$fmi$CFA, 
                                # something goes wrong with some ridge parameters
                                # for the CFA model fitting
                                out[[dt]][[i]]$fmi$semS,
                                out[[dt]][[i]]$fmi$lm)
    )
  }
  store_0 <- cbind(store_0, rowMeans(store_1))
}
ridge_range <- length(unique(out$conds$ridge))
colnames(store_0) <- names(out[[1]])
colnames(store_0) <- rep(unique(out$conds$ridge), nrow(out$conds)/ridge_range)

avg_fmi <- round(colMeans(store_0), 3)

i <- 1
j <- i+ridge_range-1
ridge_p <- NULL
for (r in 1:(length(avg_fmi)/ridge_range)) {
  ridge_p[r] <- as.numeric(names(which.min(avg_fmi[i:j])))
  i <- j+1
  j <- i+ridge_range-1
}

## Ridge that gives smallest fmi for each condition
bridge_cv <- data.frame(out$conds[!duplicated(out$conds[, -1]), -1],
                        ridge = ridge_p)
bridge_cv

## Inspect high dimensional conditions
store_0 <- NULL
for (i in 1:nrow(out$conds)) {
  store_1 <- NULL
  for (dt in 1:out$parms$dt_rep) {
    store_1 <- cbind(store_1, c(out[[dt]][[i]]$fmi$semS)
    )
  }
  store_0 <- cbind(store_0, rowMeans(store_1))
}


ridge_range <- length(unique(out$conds$ridge))
colnames(store_0) <- names(out[[1]])
colnames(store_0) <- rep(unique(out$conds$ridge), nrow(out$conds)/ridge_range)

t(store_0)
# You can see that for each paramter (columns), in the hgih dimensional conditions
# 1e-07 gives the best fmi, often 1 order of magnitude smaller than neighbour
# values.
# For low dimensional conditions. the difference is not as stark