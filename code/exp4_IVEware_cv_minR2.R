# Project:   imputeHD-comp
# Objective: Cross-validation of minR2 on EVS data
# Author:    Edoardo Costantini
# Created:   2023-03-22
# Modified:  2023-03-22
# Notes: 

# Prepare environment
rm(list = ls())
source("./init_general.R")
source("./exp4_init.R")
source("../crossvalidate/exp4_cv_minR2_init.R")

# Check where are you running IVEware is correct for OS
parms$IVEloc

# Check parameters in use
parms

# Check conditions in use
conds

# Create a cluster object:
clus <- makeCluster(5)

# Prep Environments
clusterEvalQ(cl = clus, expr = source("./init_general.R"))
clusterEvalQ(cl = clus, expr = source("./exp4_init.R"))
clusterEvalQ(cl = clus, expr = source("../crossvalidate/exp4_cv_minR2_init.R"))

# Progress report file ----------------------------------------------------

# Create a progress report file
file.create(paste0(parms$outDir, parms$report_file_name))

# Add prints
cat(
        paste0(
                "IVEWARE MINR2 CROSSVALIDATION PROGRESS REPORT",
                ## Description
                "\n",
                "Description: ", parms$description,
                "\n",
                ## Time
                "\n", "------", "\n",
                "Starts at: ", Sys.time(),
                "\n", "------", "\n"
        ),
        file = paste0(parms$outDir, parms$report_file_name),
        sep = "\n",
        append = TRUE
)

# parallel apply ---------------------------------------------------------------

# Time stamp: start
sim_start <- Sys.time()

# Run the computations in parallel on the 'clus' object:
out <- parLapply(
        cl = clus,
        X = 1:parms$dt_rep,
        fun = doRep,
        conds = conds,
        parms = parms,
        debug = FALSE,
        verbose = FALSE
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
out <- readRDS("../output/exp4_cv_IVEware_20230322_1841.rds")

# Obtain conditions with cv ridge
minR2_cv <- cvParm(
        out = out,
        cv.parm = colnames(out$conds)[1],
        exp_factors = colnames(out$conds)[2]
)

# Look at plot
minR2_cv$plot

# Look at solutions
minR2_cv$values