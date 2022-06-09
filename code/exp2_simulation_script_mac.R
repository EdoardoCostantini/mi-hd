### Title:    Imputing High Dimensional Data
### Author:   Anonymized for peer review
### Created:  2020-07-22

rm(list=ls())
source("./init_general.R")
source("./exp2_init.R")

## Data directory for storage

# Progress report file ----------------------------------------------------
file.create(paste0(parms$outDir, parms$report_file_name))

cat(paste0("SIMULATION PROGRESS REPORT",
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

  out <- mclapply(X        = 1 : parms$dt_rep,
                  FUN      = doRep,
                  conds    = conds,
                  parms    = parms,
                  debug    = FALSE,
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
  out$parms <- parms
  out$conds <- conds
  
# Save output -------------------------------------------------------------

saveRDS(out,
        paste0(parms$outDir,
               parms$results_file_name)
)