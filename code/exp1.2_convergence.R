# Project:   imputeHD-comp
# Objective: Check convergence for collinearity study
# Author:    Edoardo Costantini
# Created:   2023-03-28
# Modified:  2023-04-03
# Notes: 

# 1. Convergence for all R-based MI methods -------------------------------------

# Clean environment
rm(list = ls())

# Source scripts
source("./init_general.R")
source("./exp1.2_init.R")

# 1.2 Modify parms and conds objects -------------------------------------------

# 1.2.1 Parameters -------------------------------------------------------------

# Which imputation method to use
parms$meth_sel <- data.frame(
    DURR_la = TRUE,
    IURR_la = TRUE,
    blasso = TRUE,
    bridge = TRUE,
    MI_PCA = TRUE,
    MI_CART = TRUE,
    MI_RF = TRUE,
    stepFor = FALSE,  # <- Special method requiring a different strategy
    MI_qp = TRUE,
    MI_am = TRUE,
    MI_OP = TRUE,
    missFor = TRUE,
    mean = TRUE,
    CC = TRUE,
    GS = TRUE
)
parms$methods <- names(parms$meth_sel)[which(parms$meth_sel == TRUE)]

# What to store
parms$store <- c(
    cond = TRUE,
    dat_full = TRUE,
    dat_miss = TRUE,
    sem_EST = FALSE,
    sem_CI = FALSE,
    lm_EST = FALSE,
    lm_CI = FALSE,
    fmi = FALSE,
    miss_descrps = FALSE,
    run_time_min = FALSE,
    imp_values = TRUE # <- this is what we want for traceplots
)

# Iterations, repetitions, etc
parms$dt_rep <- 5
parms$ndt <- parms$mice_ndt <- parms$chains <- 5
parms$iters <- parms$mice_iters <- 10 # 100 iterations
parms$burnin_imp <- 0
parms$thin <- (parms$iters - parms$burnin_imp) / parms$ndt
parms$pos_dt <- (parms$burnin_imp + 1):parms$iters # candidate datasets (after convergence)
parms$keep_dt <- parms$pos_dt[seq(
    1,
    length(parms$pos_dt),
    parms$thin
)] # keep 1 dataset every thin

# For blasso
parms$chains_bl <- parms$chains
parms$iters_bl <- parms$iters
parms$burnin_imp_bl <- 0
parms$thin_bl <- (parms$iters_bl - parms$burnin_imp_bl) / parms$ndt
parms$pos_dt_bl <- (parms$burnin_imp_bl + 1):parms$iters_bl
parms$keep_dt_bl <- parms$pos_dt_bl[seq(
    1,
    length(parms$pos_dt_bl),
    parms$thin_bl
)]

# Report names
filename <- paste0(
    "exp",
    parms$exp,
    "_convergence_all_meth_",
    parms$start_time
)
parms$report_file_name <- paste0(filename, ".txt")
parms$results_file_name <- paste0(filename, ".rds")

# 1.2.2 Conditions -------------------------------------------------------------

latent <- FALSE
pm <- .3
p <- 500 # c(50, 500) # number of variables
collinearity <- c(.9)
ridge <- 0.01

# Create experimental conditions
conds <- expand.grid(
    ridge = ridge,
    p = p,
    latent = latent,
    pm = pm,
    collinearity = collinearity
)

# Bridge special parameters per condition
conds$ridge <- 0.01

# 1.3 Perform simulation -------------------------------------------------------

# Start simulation 
sim_start <- Sys.time()

# Create a cluster object:
clus <- makeCluster(parms$dt_rep)

# Export call for packages
clusterEvalQ(cl = clus, expr = source("./init_general.R"))

# Export parms and conds objects
clusterExport(
    cl = clus,
    varlist = list(
        "parms",
        "conds"
    )
)

# Progress report file
file.create(paste0(parms$outDir, parms$report_file_name))

cat(
    paste0(
        "SIMULATION PROGRESS REPORT",
        ## Description
        "\n",
        "Description: ", parms$description,
        "\n",
        ## Time
        "\n", "------", "\n",
        "Starts at: ", Sys.time(),
        "\n", "------", "\n"
    ),
    file = paste0(parms$outDir, parms$report_file_name),
    sep = "\n",
    append = TRUE
)

# Run the computations in parallel on the 'clus' object:
out <- parLapply(
    cl = clus,
    X = 1:(parms$dt_rep),
    fun = doRep,
    conds = conds,
    parms = parms
)

# Kill the cluster:
stopCluster(clus)

sim_ends <- Sys.time()

# Append Parms and conds to out object
out$parms <- parms
out$conds <- conds
out$session_info <- devtools::session_info()

# Close report file
cat(
    paste0(
        "\n", "------", "\n",
        "Ends at: ", Sys.time(), "\n",
        "Run time: ",
        round(difftime(sim_ends, sim_start, units = "hours"), 3), " h",
        "\n", "------", "\n"
    ),
    file = paste0(parms$outDir, parms$report_file_name),
    sep = "\n",
    append = TRUE
)

# Save results
saveRDS(
    out,
    paste0(
        parms$outDir,
        parms$results_file_name
    )
)

# 1.4 Analyze results ----------------------------------------------------------

# Load data
out_cnv <- readRDS("../output/exp1_2_convergence_all_meth_20230403_1027.rds") # colli .9

# What to show
iters_range <- 1:100 # which set of iterations
y_range <- c(3, 6)
exp_dat <- 2 # which data replication (10 possibilities)

# Check run is correct for all involved methods
data.frame(
    Chians = c(
        all = out_cnv$parms$chains, 
        blasso = out_cnv$parms$chains_bl,
        MICE = out_cnv$parms$mice_ndt
    ),
    Iters = c(
        all = out_cnv$parms$iters,
        blasso = out_cnv$parms$iters_bl,
        MICE = out_cnv$parms$mice_iters
    )
)

# DURR
mean_traceplot(
    out = out_cnv,
    method = "DURR_la",
    dat = exp_dat,
    y_range = y_range,
    iters = iters_range
)

# IURR
mean_traceplot(
    out = out_cnv,
    method = "IURR_la",
    dat = exp_dat,
    y_range = y_range,
    iters = iters_range
)

# blasso
mean_traceplot(
    out = out_cnv,
    method = "blasso",
    dat = exp_dat,
    y_range = c(3, 7), 
    iters = iters_range
)

# bridge
mean_traceplot(
    out = out_cnv,
    method = "bridge",
    dat = exp_dat,
    y_range = c(-5, 15),
    iters = iters_range
)

# MI_PCA
mean_traceplot(
    out = out_cnv,
    method = "MI_PCA",
    dat = exp_dat,
    y_range = y_range, iters = iters_range
)

# MI_CART
mean_traceplot(
    out = out_cnv,
    method = "MI_CART",
    dat = exp_dat,
    y_range = y_range, iters = iters_range
)

# MI_RF
mean_traceplot(
    out = out_cnv,
    method = "MI_RF",
    dat = exp_dat,
    y_range = y_range, iters = iters_range
)

# MI_qp
mean_traceplot(
    out = out_cnv,
    method = "MI_qp",
    dat = exp_dat,
    y_range = y_range,
    iters = iters_range
)

# MI_am
mean_traceplot(out_cnv,
    method = "MI_am",
    dat = exp_dat,
    y_range = y_range,
    iters = iters_range
)

# MI_OP
mean_traceplot(
    out = out_cnv,
    method = "MI_OP",
    dat = exp_dat,
    y_range = y_range,
    iters = iters_range
)

# 2. Convergence for IVEware ---------------------------------------------------
# TODO: this is not ready for simulation
# Clean environment
rm(list = ls())

# Source scripts
source("./init_general.R")
source("./exp1.2_init.R")

# 2.2 Modify parms and conds objects -------------------------------------------

# 2.2.1 Parameters -------------------------------------------------------------

# Which imputation method to use
parms$meth_sel <- data.frame(
    DURR_la = FALSE,
    IURR_la = FALSE,
    blasso = FALSE,
    bridge = FALSE,
    MI_PCA = FALSE,
    MI_CART = FALSE,
    MI_RF = FALSE,
    stepFor = TRUE, # <- method of interest
    MI_qp = FALSE,
    MI_am = FALSE,
    MI_OP = FALSE,
    missFor = TRUE,
    mean = TRUE,
    CC = TRUE,
    GS = TRUE
)
parms$methods <- names(parms$meth_sel)[which(parms$meth_sel == TRUE)]

# What to store
parms$store <- c(
    cond = TRUE,
    dat_full = TRUE,
    dat_miss = TRUE,
    sem_EST = FALSE,
    sem_CI = FALSE,
    lm_EST = FALSE,
    lm_CI = FALSE,
    fmi = FALSE, 
    miss_descrps = FALSE,
    run_time_min = FALSE,
    imp_values = TRUE #  <- this is what we want for convergence checks
)

# Iterations, repetitions, etc
parms$dt_rep <- 20

# Report names
filename <- paste0(
    "exp",
    parms$exp, "_",
    "cv_IVEware_",
    parms$start_time
)
parms$report_file_name <- paste0(filename, ".txt")
parms$results_file_name <- paste0(filename, ".rds")

# 2.2.2 Conds ------------------------------------------------------------------

p <- 500
pm <- .3
iters <- c(1, 5, 10, 20, 40, 80, 160, 240, 320)
seed <- c(2023, 3202, 1610, 4004, 9876)
conds <- expand.grid(p = p, pm = pm, iters = iters, seed = seed)

# 2.3 Perform simulation -------------------------------------------------------

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
        minR2 = 0.01,
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

# 2.4 Analyze results ----------------------------------------------------------

# Read results
conv.out <- readRDS("../output/")

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