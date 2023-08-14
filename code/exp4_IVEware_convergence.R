# Project:   imputeHD-comp
# Objective: Check the convergence of IVEware through density plots on EVS set up
# Author:    Edoardo Costantini
# Created:   2023-03-28
# Modified:  2023-03-28
# Notes:

# Prepare environment
rm(list = ls())
source("./init_general.R")
source("./exp4_init.R")
source("../convergence/exp4_ccheck_init_IVEware.R")

# Set a seed
set.seed(1234)

# Source the data
data_source <- readRDS("../data/exp4_EVS2017_full.rds")$full

# Generated a bootstrap sample
Xy <- data_source[
    sample(1:nrow(data_source),
        conds[1, ]$n,
        replace = TRUE
    ),
]

# Impose missing values
Xy_mis <- imposeMiss_evs(Xy, parms, conds[1, ])

# Missing data
miss_descrps <- colMeans(is.na(Xy_mis)[, parms$z_m_id])

# Starting time stamp
when <- format(Sys.time(), "%Y%m%d_%H%M")

# Create an empty object to store the results
store_imps <- list()

# Define condition
for (i in 1:nrow(conds)) {
    # Which condition
    cond <- conds[i, ]

    # Set condition seed
    set.seed(cond$seed)

    # Perform imputation
    store_imps[[i]] <- impute_IVEware(
        Z = Xy_mis,
        minR2 = 0.001,
        rep_status = i,
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
conv.out <- readRDS("../output/exp4_cv_IVEware_20230328_1544.rds")

# For every imputed variable

pdf(file = "../output/graphs/exp4_cv_IVEware_20230328_1544.pdf", width = 20, height = 20)
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
            main = paste0("Density plot for observed values and imputed imputed on ", variable),
            xlab = "",
            ylim = c(0, max(maxDens))
        )

        # Add densities of missing values
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

dev.off()
