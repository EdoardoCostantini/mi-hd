### Title:    Imputing High Dimensional Data
### Author:   Anonymized for peer review
### Created:  2020-05-19

rm(list=ls())
source("./init_general.R")
source("./exp4_init.R")

## Create a cluster object:
clus <- makeCluster(5) # 30 on blade

## Two different ways to source a script on the worker nodes:
clusterEvalQ(cl = clus, expr = source("./init_general.R"))
clusterEvalQ(cl = clus, expr = source("./exp4_init.R"))

# Progress report file ----------------------------------------------------
file.create(paste0(parms$outDir, parms$report_file_name))

cat(paste0("SIMULATION PROGRESS REPORT",
           ## Description
           "\n",
           "Description: ", parms$description, 
           "\n",
           ## Time
           "\n", "------", "\n",
           "Starts at: ", Sys.time(),
           "\n", "------", "\n" ),
    file = paste0(parms$outDir, parms$report_file_name),
    sep = "\n",
    append = TRUE)

# mcApply parallel --------------------------------------------------------

sim_start <- Sys.time()

## Run the computations in parallel on the 'clus' object:
out <- parLapply(cl = clus, 
                 X = 1 : parms$dt_rep,
                 fun = doRep, 
                 conds = conds, 
                 parms = parms,
                 debug = FALSE,
                 verbose = FALSE)

## Kill the cluster:
stopCluster(clus)

sim_ends <- Sys.time()

# Attach Info Objects
out$parms <- parms
out$conds <- conds
out$session_info <- session_info()

## Close report file
cat(paste0("\n", "------", "\n",
           "Ends at: ", Sys.time(), "\n",
           "Run time: ",
           round(difftime(sim_ends, sim_start, units = "hours"), 3), " h",
           "\n", "------", "\n"),
    file = paste0(parms$outDir, parms$report_file_name),
    sep = "\n",
    append = TRUE)

# Save output -------------------------------------------------------------

saveRDS(out,
        paste0(parms$outDir,
               parms$results_file_name)
)
