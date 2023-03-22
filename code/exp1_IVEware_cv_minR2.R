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

# Create a progress report file
file.create(paste0(parms$outDir, parms$report_file_name))
cat(
        paste0(
                "IVEWARE MINR2 CROSSVALIDATION PROGRESS REPORT",
                "\n",
                "Description: ", parms$description, "\n",
                "\n", "------", "\n",
                "Starts at: ", Sys.time(),
                "\n", "------", "\n"
        ),
        file = paste0(parms$outDir, parms$report_file_name),
        sep = "\n",
        append = TRUE
)

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

# Take note of end in progress report
cat(
        paste0(
                "\n", "------", "\n",
                "Ends at: ", Sys.time(), "\n",
                "Run time: ",
                round(difftime(sim_ends, sim_start, units = "hours"), 3), " h",
                "\n", "------", "\n"
        ),
        file = paste0(parms$outDir, parms$report_file_name),
        sep = "\n",
        append = TRUE
)

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