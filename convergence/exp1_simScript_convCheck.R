### Title:    Convergence Checks: Initialization Script for Exp 1 
### Project:  Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-07-03
### Notes:    Runs the simulation for convergence checks (few repetitions, 
###           more chains)
### IMPORTANT: Working directory: ~/imputeHD-comp/code

rm(list = ls())
source("./init_general.R")
source("../convergence/exp1_init_convCheck.R")

## Data directory for storage

# Progress report file ----------------------------------------------------
file.create(paste0(parms$outDir, parms$report_file_name))

cat(paste0("CONVERGENCE CHECK PROGRESS REPORT",
           "\n",
           "Description: ", parms$description, "\n",
           "\n", "------", "\n",
           "Starts at: ", Sys.time(),
           "\n", "------", "\n" ),
    file = paste0(parms$outDir, parms$report_file_name),
    sep = "\n",
    append = TRUE)

# mcApply parallel --------------------------------------------------------

sim_start <- Sys.time()

out_cnv <- mclapply(X        = 1 : parms$dt_rep,
                FUN      = doRep,
                conds    = conds,
                parms    = parms,
                mc.cores = ( 10 ) )

sim_ends <- Sys.time()

cat(paste0("\n", "------", "\n",
           "Ends at: ", Sys.time(), "\n",
           "Run time: ", 
           round(difftime(sim_ends, sim_start, units = "hours"), 3), " h",
           "\n", "------", "\n"),
    file = paste0(parms$outDir, parms$report_file_name),
    sep = "\n",
    append = TRUE)

# Attach parm object
out_cnv$parms <- parms

# Save output -------------------------------------------------------------

saveRDS(out_cnv,
        paste0(parms$outDir,
               parms$results_file_name)
)