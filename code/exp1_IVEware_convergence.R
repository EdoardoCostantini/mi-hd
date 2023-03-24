# Project:   imputeHD-comp
# Objective: Check the convergence of IVEware through density plots
# Author:    Edoardo Costantini
# Created:   2023-03-24
# Modified:  2023-03-24
# Notes: 

# Prepare environment
rm(list = ls())
source("./init_general.R")
source("./exp1_init.R")
source("../convergence/exp1_ccheck_init_IVEware.R")

# Set a seed
set.seed(1234)

# Create an empty object to store the results
store_imps <- list()

# Define condition
for (i in 1:nrow(conds)) {
    # Which condition
    cond <- conds[i, ]

    # Generated data
    Xy <- simData_exp1(cond, parms)

    # Impose missing values
    Xy_mis <- imposeMiss(Xy, parms, cond)
    Xy_mis <- cbind(
        Xy_mis[, parms$z_m_id],
        Xy_mis[, -which(colnames(Xy_mis) %in% parms$z_m_id)]
    )

    # Define matrix index of observed values
    O <- !is.na(Xy_mis)

    # Define the number of iterations to check
    parms$iters <- cond$iters

    # Define the number of multiply imputed datasets
    parms$ndt <- 30

    # Perform imputation
    store_imps[[i]] <- impute_IVEware(
        Z = Xy_mis,
        minR2 = 0.001,
        rep_status = 1,
        iters = cond$iters,
        perform = TRUE,
        parms = parms
    )
}

# Save results
saveRDS(
    store_imps,
    paste0(
        "exp1_conv_IVEware_",
        format(Sys.time(), "%Y%m%d_%H%M"),
        ".rds"
    )
)

# Read results
store_imps <- readRDS("")

# For every imputed variable
par(mfrow = c(3, 2))

# Loop over conditions
for (condition in 1:nrow(conds)) {
    # Look over variables
    for (variable in parms$z_m_id) {
        # Plot baseline
        plot(
            density(na.omit(Xy_mis[, variable])),
            lwd = 2,
            main = paste0("Density plot for observed and imputed ", variable),
            xlab = ""
        )
        # Add densities for
        lapply(1:parms$ndt, function(j) {
            lines(
                density(store_imps[[condition]]$dats[[j]][, variable]),
                col = "blue",
                lwd = 1
            )
        })
        # Add a shared title
        mtext(
            paste0("Iterations = ", conds$iters[condition]),
            side = 1, # bottom placement
            line = -2,
            outer = TRUE
        )
    }
}
