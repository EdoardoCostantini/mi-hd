# Object: Stores infomration on the session that will be run on Lisa
# Notes:  This script should be called by a prepartion script that you run 
#         on Lisa

# Prepare Environment
rm(list = ls())
source("./init_general.R")
source("./exp5_init.R")

# Output Directory from Terminal inputs
args <- commandArgs(trailingOnly = TRUE)
parms$outDir <- args[1]

# Create Empty storing object
out <- list(parms = parms,
            conds = conds,
            session_info = session_info())

# Save it in the root
saveRDS(out,
        paste0(parms$outDir,
               parms$runInfo)
)
