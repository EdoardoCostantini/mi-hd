### Title:    Convergence Checks: Single Imputation for pre-processing for addendum
### Project:  Imputing High Dimensional Data Addendum
### Author:   Edoardo Costantini
### Created:  2021-01-29
### Notes:    Generates data and check what values should be used
###           for single imputation in the preprocessing phase
### IMPORTANT: Working directory: ~/imputeHD-comp/code

rm(list=ls())
source("./init_general.R")
source("./exp5_init.R")

# Single Run function
runCell <- function(rp, cond) {
    ## Function for parallel runs of simulation > impose miss > single imp
    ## The results are then displayed with a different function
    ## Example Inputs
    # cond = conds[1, ]
    
    ## Body
    simData_list <- simData_lv(parms, cond = cond)
    Xy     <- simData_list$dat
    Xy_mis <- imposeMiss_lv_MD(simData_list, 
                               parms, 
                               cond, 
                               plot = TRUE)
    
    O <- !is.na(Xy_mis) # matrix index of observed values
    
    # Prep Data Perform Single Imputation
    prep_Si_out <- prep_SI(dt_in = Xy_mis,
                           model.var = parms$z_m_id,
                           # for internal mice run
                           m = 5,
                           maxit = 200,
                           ridge = cond$ridge)
    
    # Outcome
    return(prep_Si_out$mids_obj)
}

# Example Run of a single Check
# output <- runCell(rp = 1, cond = conds[1, ])
# plot(output,
#      y = colnames(output$data)[10:19],
#      layout = c(2, 10),
#      type = "l")

# Actual Check
# I want to see for a sample of vairables to be imputed with single 
# imputation, when more or less convergence is reached.
# To do so, we can study MI convergence on 10 different replication
# of the data generation procedure of interest and then compare how
# the traceplots change across these repetitions

## Create a cluster object:
clus <- makeCluster(10)

## Source scripts on the worker nodes
clusterEvalQ(cl = clus, expr = source("./init_general.R"))
clusterEvalQ(cl = clus, expr = source("./exp5_init.R"))

## Run the computations in parallel on the 'clus' object:
sim_start <- Sys.time()

out <- parLapply(cl = clus, 
                 X = 1 : 10,
                 fun = runCell, 
                 cond = conds[nrow(conds), ])

sim_end <- Sys.time()
sim_end - sim_start

## Kill the cluster:
stopCluster(clus)

## Check plots latent variable by latent variable
var_lists <- names(which(rowSums(out[[1]]$predictorMatrix) != 0))

# Look at the hand written table in your notion notebook w/ heading 
# "Defining the Set up: Simulation or Resampling study"
# to see why you chose these variables
item_list <- list(LV4 = (4-1) * 5+(1:5), # Block B (no-MAR high cor latent var)
                  LV5 = (5-1) * 5+(1:5), # Block A (no-MAR low corr latent var)
                  LV6 = (6-1) * 5+(1:5), # Blcok B (no-MAR low corr latent var)
                  LV7 = (7-1) * 5+(1:5), # Blcok C no-MAR low corr latent var
                  LV8 = (8-1) * 5+(1:5), # Blcok C no-MAR low corr latent var
                  LV70 = (70-1) * 5+(1:5),  # Junk no-MAR low corr latent var
                  LV100 = (100-1) * 5+(1:5)  # Junk no-MAR low corr latent var
                  )

pdf(file = "../output/exp5_ccheck_Si_20210326_1750.pdf",
    width = 15, height = 15)
lapply(seq_along(item_list), function(i){
    # i <- 2
    lapply(seq_along(out), function(x){
        # x <- 2
        plot(out[[x]],
             main = paste0(names(item_list)[i],
                           " Rep", x),
             y = colnames(out[[x]]$data)[item_list[[i]]],
             layout = c(2, 5),
             type = "l")
    })
})
dev.off()

## Save results
saveRDS(out, paste0("parms$outDir", 
                    "exp5_ccheck_S_",
                    format(Sys.time(), "%Y%m%d_%H%M"),
                    ".rds"))

## Read Old Results
out <- readRDS("../output/exp5_ccheck_Si_20210326_1750.rds")
length(out)
