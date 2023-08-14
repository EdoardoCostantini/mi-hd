# Project:   imputeHD-comp
# Objective: Check the convergence of IVEware through density plots
# Author:    Edoardo Costantini
# Created:   2023-03-24
# Modified:  2023-03-27
# Notes: 

# Prepare environment
rm(list = ls())
source("./init_general.R")
source("./exp1_init.R")
source("../convergence/exp1_ccheck_init_IVEware.R")

# Set a seed
set.seed(1234)

# Generated data
Xy <- simData_exp1(conds[1, ], parms)

# Impose missing values
Xy_mis <- imposeMiss(Xy, parms, conds[1, ])
Xy_mis <- cbind(
    Xy_mis[, parms$z_m_id],
    Xy_mis[, -which(colnames(Xy_mis) %in% parms$z_m_id)]
)

# Starting time stamp
when <- format(Sys.time(), "%Y%m%d_%H%M")

# Create an empty object to store the results
store_imps <- list()

# Define condition
for (i in 1:nrow(conds)) {
    # Set condition seed
    set.seed(conds[i, "seed"])

    # Which condition
    cond <- conds[i, ]

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
    list(
        original.data = Xy_mis,
        imputed.data = store_imps,
        parms = parms,
        conds = conds
    ),
    paste0(
        "../output/", parms$results_file_name
    )
)

# Read results
conv.out <- readRDS("../output/exp1_conv_IVEware_20230327_1143.rds")

# For every imputed variable
pdf(file = "../output/graphs/exp1_conv_IVEware_20230327_1143.pdf", width = 20, height = 20)
par(mfrow = c(6, 2))

# Loop over conditions
for (condition in 1:nrow(conv.out$conds)) {
    # Look over variables
    for (variable in conv.out$parms$z_m_id) {

        # Identify imputed values
        imps <- is.na(conv.out$original.data[, variable])

        # Find the highest point for ylim
        maxDens <- sapply(1:conv.out$parms$ndt, function(j) {
            object <- density(conv.out$imputed.data[[1]]$dats[[j]][imps, variable])
            max(object$y)
        })

        # Plot baseline
        plot(
            density(na.omit(conv.out$original.data[, variable])),
            lwd = 2,
            main = paste0("Density plot for observed and imputed ", variable),
            xlab = ""
        )
        # Add densities for
        lapply(1:conv.out$parms$ndt, function(j) {
            lines(
                density(conv.out$imputed.data[[condition]]$dats[[j]][imps, variable]),
                col = "blue",
                lwd = 1
            )
        })

        # Plot baseline
        plot(
            density(na.omit(conv.out$original.data[, variable])),
            lwd = 2,
            main = paste0("Density plot for observed and imputed ", variable),
            xlab = "",
            ylim = c(0, max(maxDens))
        )

        # Add densities of imputed variable
        lapply(1:conv.out$parms$ndt, function(j) {
            lines(
                density(conv.out$imputed.data[[condition]]$dats[[j]][, variable]),
                col = "red",
                lwd = 1
            )
        })

        # Add a shared title
        mtext(
            paste0(
                "Seed = ", conv.out$conds$seed[condition], "; ",
                "Iterations = ", conv.out$conds$iters[condition]
                ),
            side = 1, # bottom placement
            line = -2,
            outer = TRUE
        )
    }
}

# Close pdf storing
dev.off()
