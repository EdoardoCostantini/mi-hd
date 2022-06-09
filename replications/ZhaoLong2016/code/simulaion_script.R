### Title:    Replication Zhao Long 2016 - Simulation script
### Author:   Anonymized for peer review
### Created:  2020-03-03
### Modified: 2020-03-03

rm(list=ls())
source("./init.R")

## Data directory for storage

# Progress report file ----------------------------------------------------
file.create(paste0(parms$outDir, parms$report_file_name))

cat(paste0("SIMULATION PROGRESS REPORT",
           "\n",
           "Date: ", Sys.Date(),
           "\n",
           "Description: ", parms$description,
           "\n"),
    file = paste0(parms$outDir, parms$report_file_name),
    sep = "\n",
    append = TRUE)

# mcApply parallel --------------------------------------------------------

out <- mclapply(X        = 1 : parms$dt_rep,
                FUN      = doRep,
                conds    = conds,
                parms    = parms,
                mc.cores = ( 10 ) )

# Save output -------------------------------------------------------------

saveRDS(out,
        paste0(parms$outDir,
               parms$results_file_name)
)
