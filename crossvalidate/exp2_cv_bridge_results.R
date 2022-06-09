### Title:    Conergence Check Blasso
### Project:  Imputing High Dimensional Data
### Author:   Anonymized for peer review
### Created:  2020-07-03
### Notes:    The goal is to find out how many iterations should be used
###           in the full study for each method.

rm(list = ls())
source("./init_general.R")

out <- readRDS("../output/exp2_cv_bridge_20200812_1132.rds")

# Extract average fmi for a given model type ------------------------------
conds_bridge <- bridge_cv(out, 
                          mods = names(out[[1]][[1]]$fmi)[-2] # CFA failed fitting
                          )
conds_bridge

# Inspect high dimensional conditions -------------------------------------

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