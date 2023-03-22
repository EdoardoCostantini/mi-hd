# Title:    Simulation script for ridge paramter cross-validation
# Project:  Imputing High Dimensional Data
# Author:   Edoardo Costantini
# Created:  2020-08-25
# Modified: 2023-03-22

# Prepare environment
rm(list = ls())
source("./init_general.R")
source("./exp1_init.R")

# Check where are you running IVEware is correct for OS
parms$IVEloc

# Time stamp: start
sim_start <- Sys.time()

# Create a cluster object:
clus <- makeCluster(5)

# Prep Environments
clusterEvalQ(cl = clus, expr = source("./init_general.R"))
clusterEvalQ(cl = clus, expr = source("./exp1_init.R"))
clusterEvalQ(cl = clus, expr = source("../crossvalidate/exp1_cv_minR2_init.R"))

# Run the computations in parallel on the 'clus' object:
out <- parLapply(
        cl = clus,
        X = 1:parms$dt_rep,
        fun = doRep,
        conds = conds,
        parms = parms
)

# Kill the cluster:
stopCluster(clus)

# Time stamp: end
sim_ends <- Sys.time()

# Append parms and conds to out object
out$parms <- parms
out$conds <- conds
out$session_info <- devtools::session_info()

# Save results
saveRDS(
        out,
        paste0(
                parms$outDir,
                parms$results_file_name
        )
)

# Read the results again
out <- readRDS("../output/exp1_cv_IVEware_20230321_1748.rds")

# Obtain conditions with cv ridge
minR2_cv <- cvParm(
        out = out,
        cv.parm = colnames(out$conds)[1],
        exp_factors = colnames(out$conds)[c(2, 4)]
)

# Look at plot
minR2_cv$plot

# Look at solutions
minR2_cv$values