### Title:    Run imputations (parallel)
### Author:   Anonymized for peer review
### Created:  2020-02-13
### Modified: 2020-02-15
### Notes:    this version using mcapply allows you

library(parallel) # detectCores(); makeCluster()

source("./init_mcapply.R")

## Data directory for storage

outDir <- "../output/"

# mcApply parallel --------------------------------------------------------

out <- mclapply(X        = 1 : parms$dt_rep,
                FUN      = singleRun,
                chains   = parms$chains,
                iter     = parms$iters, 
                parms    = parms,
                mc.cores = ( detectCores()-1 ) )

# Save output -------------------------------------------------------------

saveRDS(out,
        paste0(outDir,
               "pooled_BR2010-mc-",
               format(Sys.time(), "%Y%m%d_%H%M"),
               ".rds")
)