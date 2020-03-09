### Title:    Replication Zhao Long 2016 - Simulation script
### Author:   Edoardo Costantini
### Created:  2020-03-03
### Modified: 2020-03-03

library(parallel) # detectCores(); makeCluster()

source("./init.R")

## Data directory for storage
outDir <- "../output/"

# mcApply parallel --------------------------------------------------------

out <- mclapply(X        = 1 : parms$dt_rep,
                FUN      = doRep,
                conds    = conds,
                parms    = parms,
                mc.cores = ( 10 ) )

# Save output -------------------------------------------------------------

saveRDS(out,
        paste0(outDir,
               "pooled_ZL2016-mc-",
               format(Sys.time(), "%Y%m%d_%H%M"),
               ".rds")
)
