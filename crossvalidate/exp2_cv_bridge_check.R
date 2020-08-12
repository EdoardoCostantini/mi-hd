### Title:    Conergence Check Blasso
### Project:  Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-07-03
### Notes:    The goal is to find out how many iterations should be used
###           in the full study for each method.

rm(list = ls())
source("./init_general.R")
source("./exp2_init.R.R")
source("../crossvalidate/exp2_cv_bridge_init.R")

sim_start <- Sys.time()

# Create a cluster object:
clus <- makeCluster(parms$dt_rep)

# Prep Environments
clusterEvalQ(cl = clus, expr = source("./init_general.R"))
clusterEvalQ(cl = clus, expr = source("./init_exp2.R"))
clusterEvalQ(cl = clus, expr = source("../crossvalidate/cv_bridge_init_exp2.R"))

## Run the computations in parallel on the 'clus' object:
out <- parLapply(cl = clus, 
                 X = 1 : parms$dt_rep,
                 fun = doRep, 
                 conds = conds, 
                 parms = parms)

## Kill the cluster:
stopCluster(clus)

sim_ends <- Sys.time()

out$parms <- parms
out$conds <- conds

# Extract average fmi for a given model type ------------------------------
store_0 <- NULL
for (i in 1:nrow(out$conds)) {
  store_1 <- NULL
  for (dt in 1:out$parms$dt_rep) {
    store_1 <- cbind(store_1, c(out[[dt]][[i]]$fmi$semR,
                                out[[dt]][[i]]$fmi$CFA,
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

## Show result
bridge_cv <- data.frame(out$conds[!duplicated(out$conds[, -1]), -1],
                        ridge = ridge_p)


saveRDS(bridge_cv,
        paste0(parms$outDir,
               parms$results_file_name)
)
