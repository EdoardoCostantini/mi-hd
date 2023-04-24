# Project:   imputeHD-comp
# Objective: Check different aspect of the code base
# Author:    Edoardo Costantini
# Created:   2020-05-19
# Modified:  2023-04-21
# Notes: 

rm(list=ls())
source("./init_general.R")
source("./exp2_init.R")

# SEM Model convergence ---------------------------------------------------

cond <- data.frame(p = 500, pm = .3)
# fit_gs <- list()
# fit_cc <- list()

set.seed(20200625)

for (i in 1:5) {
  print(i)
  Xy <- genDt_mvn(cond, parms)
  Xy_mis <- imposeMiss(Xy, parms, cond)
  Xy_mis <- cbind(Xy_mis[parms$z_m_id], Xy_mis[-parms$z_m_id])
  
  O <- !is.na(Xy_mis) # matrix index of observed values
  miss_descrps <- colMeans(!O[, 1:parms$zm_n]) 
  
  sem_cnvg <- FALSE
  
  while (sem_cnvg == FALSE) {
    imp_PCA <- impute_PCA(Z = Xy_mis, parms = parms)
    imp_PCA_fit <- tryCatch({fit_sat_model(imp_PCA$dats)}, 
                          error = function(report) {
                            err <- paste0("Original Error: ", report)
                            return(err)
                          },
                          warning = function(report) {
                            err <- paste0("Original Warning: ", report)
                            return(err)
                          })
    sem_cnvg <- !is.character(imp_PCA_fit)
  }
  
  # For each imp method, pool estimates across the m datasets
  get_pool_EST(imp_PCA_fit)
  get_pool_CI(imp_PCA_fit)
}

which(sapply(fit_cc, class) == "try-error")

sapply(fit_cc, get_pool_EST)


summa_models <- lapply(X = fit_gs,
                       FUN = function(x) parameterEstimates(x))



# Linear Model Estiamtes --------------------------------------------------

rm(list=ls())
source("./init_general.R")
source("./exp2_init.R")

reps <- 250
cond <- conds[1,]
pm_store <- matrix(NA, nrow = reps, ncol = (length(parms$lm_model)-1))


for (r in 1:reps) {
  dt_full <- simData_exp1(cond, parms)
  lm <- lm(z1 ~ -1 + z2 + z3 + z4 + z7 + z8, as.data.frame(scale(dt_full)))
  pm_store[r,] <- lm$coefficients
  print(paste0(round(r/reps*100, 1), "% done"))
}

round(colMeans(pm_store), 3)

# Replicability -----------------------------------------------------------

# Run cells with seeds

rm(list=ls())
source("./init.R")

rp <- 1

set.seed(1234)
runCell_res1 <- runCell(cond = conds[1, ], 
                        parms = parms,
                        rep_status = rp)
set.seed(1234)
runCell_res2 <- runCell(cond = conds[1, ], 
                        parms = parms,
                        rep_status = rp)

all.equal(runCell_res1, runCell_res2)

# mcapply

rm(list=ls())
source("./init.R")

out_1 <- mclapply(X        = 1 : parms$dt_rep,
                  FUN      = doRep,
                  conds    = conds[1,],
                  parms    = parms,
                  mc.cores = ( 10 ) )

out_2 <- mclapply(X        = 1 : parms$dt_rep,
                  FUN      = doRep,
                  conds    = conds[1,],
                  parms    = parms,
                  mc.cores = ( 10 ) )

all.equal(out_1[[1]]$cond_50_0.1$all_EST, out_2[[1]]$cond_50_0.1$all_EST)
all.equal(out_1[[2]]$cond_50_0.1$all_CI, out_2[[2]]$cond_50_0.1$all_CI)
all.equal(out_1[[2]]$cond_50_0.1$imp_values, out_2[[2]]$cond_50_0.1$imp_values)

# parapply

rm(list=ls())
source("./init.R")

clus <- makeCluster(10)
clusterEvalQ(cl = clus, expr = source("./init.R"))

parapply1 <- parLapply(cl = clus, 
                 X = 1 : parms$dt_rep,
                 fun = doRep, 
                 conds = conds, 
                 parms = parms)

stopCluster(clus)

clus <- makeCluster(10)
clusterEvalQ(cl = clus, expr = source("./init.R"))

parapply2 <- parLapply(cl = clus, 
                       X = 1 : parms$dt_rep,
                       fun = doRep, 
                       conds = conds, 
                       parms = parms)

stopCluster(clus)

all.equal(parapply1[[1]]$cond_50_0.1$all_EST, parapply2[[1]]$cond_50_0.1$all_EST)
all.equal(parapply1[[2]]$cond_50_0.1$all_CI, parapply2[[2]]$cond_50_0.1$all_CI)
all.equal(parapply1[[2]]$cond_50_0.1$imp_values, parapply2[[2]]$cond_50_0.1$imp_values)
all.equal(parapply1[[5]]$cond_50_0.1$imp_values, parapply2[[5]]$cond_50_0.1$imp_values)
all.equal(parapply1[[6]]$cond_50_0.1$imp_values, parapply2[[6]]$cond_50_0.1$imp_values)
all.equal(parapply1[[8]]$cond_50_0.1$imp_values, parapply2[[8]]$cond_50_0.1$imp_values)

parapply1[[8]]

# Effects of missingness --------------------------------------------------
# Satuerated Model Estimation of Means, Variances and Covariances
# together with the CC (listwise deletion) analysis

rm(list=ls())
source("./init_general.R")
source("./exp2_init.R")
set.seed(1234)

cond <- conds[3,]

# Repeat data generation and estimations
reps <- 1e3
mu_fl <- mu_ms <- var_fl <- var_ms <- 
  matrix(NA, nrow = reps, ncol = length(parms$z_m_id))
cov_fl <- cov_ms <- 
  matrix(NA, nrow = reps, ncol = ncol(combn(parms$z_m_id, 2)))

for (r in 1:reps) {
  print(r/reps*100)
  
  Xy <- simData_exp1(cond, parms)
  Xy_mis  <- imposeMiss(Xy, parms, cond)
  O <- !is.na(Xy_mis) # matrix index of observed values
  
  sem_sndt <- lapply(list(GS      = Xy,
                          CC      = Xy_mis), # default is listwise deletion
                     sem, model = parms$lav_model, likelihood = "wishart")
  sem_par <- sem_EST(sem_sndt)
  
  # Means
  mu_fl[r, ] <- sem_par[1:parms$zm_n, "GS"]
  mu_ms[r, ] <- sem_par[1:parms$zm_n, "CC"]
  
  # Variances
  var_fl[r, ] <- sem_par[(parms$zm_n+1):(parms$zm_n*2), "GS"]
  var_ms[r, ] <- sem_par[(parms$zm_n+1):(parms$zm_n*2), "CC"]
  
  # Covariances
  cov_fl[r, ] <- sem_par[-(1:(parms$zm_n*2)), "GS"]
  cov_ms[r, ] <- sem_par[-(1:(parms$zm_n*2)), "CC"]
}

# Effect of missingness on analysis

# MCMC Estimates
  out_mu <- data.frame(full = round( colMeans(mu_fl), 3),
                       miss = round( colMeans(mu_ms), 3))
  out_va <- data.frame(full = round( colMeans(var_fl), 3),
                       miss = round( colMeans(var_ms), 3))
  out_cv <- data.frame(full = round( colMeans(cov_fl), 3),
                       miss = round( colMeans(cov_ms), 3))

# Bias (in terms of percentage of true value)
  BPR_mu <- round(abs(out_mu$full - out_mu$miss)/out_mu$full*100, 0)
  BPR_va <- round(abs(out_va$full - out_va$miss)/out_va$full*100, 0)
  BPR_cv <- round(abs(out_cv$full - out_cv$miss)/out_cv$full*100, 0)
  
# Obtain names to show things better
  fit <- sem(Xy, model = parms$lav_model, likelihood = "wishart")
  nm_par <- apply(parameterestimates(fit)[,1:3], 1, paste0, collapse = "")
  results <- data.frame(BPR=c(BPR_mu,BPR_va,BPR_cv))
  rownames(results) <- nm_par
  results
  
# Effects of missingness --------------------------------------------------

  rm(list=ls())
  source("./init_general.R")
  source("./exp2_init.R")
  parms$n <- 1e3
  
  # Gen Data
  set.seed(20200805)
  cond <- conds[3, ]
  Xy <- simData_exp1(cond, parms)
  Xy_mis <- imposeMiss(Xy, parms, cond)
  Xy_mis <- cbind(Xy_mis[, parms$z_m_id], Xy_mis[, -parms$z_m_id])
  
  # Visualize
  MI <- mice(data = Xy_mis[, 1:20], m = 50, maxit = 10,
             method = "norm")
  MI_dt <- complete(MI, "all")
  # Plots
  data.frame(i = colnames(Xy)[1:20], 
             j = colnames(Xy_mis)[1:20],
             imp = colnames(MI_dt$`1`))
  
  i <- 1; j <- 1
  
  plot(density(Xy[, i]), col = "black",
       main = "Density",
       xlab = "Z1", ylab = "",
       xlim = c(0, 10), ylim = c(0, .25), lwd = 3)
  # MAR: MI fixes the missingness
  lines(density(Xy_mis[, j], na.rm = TRUE), col = "darkorange", lwd = 3)
  lapply(MI_dt, function(x) lines(density(x[, j]), 
                                  lty = 2,
                                  col = "darkorange")
  )
  
# Ridge Penalty choice ----------------------------------------------------

  rm(list=ls())
  source("./init_general.R")
  source("./exp2_init.R")
  
  # Save a condition of low and high dimensionality from experiment 1
  cond_ld <- conds[1, ]
  cond_hd <- conds[4, ]
  
  # 1. Generate the data
  Xy_ld <- simData_exp1(cond_ld, parms)
  Xy_hd <- simData_exp1(cond_hd, parms)
  
  # 2. Impose miss
  Xy_mis_ld <- imposeMiss(Xy_ld, parms, cond_ld)
  Xy_mis_ld <- cbind(Xy_mis_ld[, parms$z_m_id], Xy_mis_ld[, -parms$z_m_id])
  
  Xy_mis_hd <- imposeMiss(Xy_hd, parms, cond_hd)
  Xy_mis_hd <- cbind(Xy_mis_hd[, parms$z_m_id], Xy_mis_hd[, -parms$z_m_id])
  
  # 4. Perform imputation w/ w/o ridge penalty based on chunk of impude_BRIDGE
  Zm_ld <- init_dt_i(Xy_mis_ld, missing_type(Xy_mis_ld)) # initialize data for each chain
  Zm_hd <- init_dt_i(Xy_mis_hd, missing_type(Xy_mis_hd))
  
  # Obtain posterior draws for all paramters of interest
  ridge_trial <- 1e-5
  j <- 1
  
  pdraw_ld <- .norm.draw(y       = Zm_ld[, j],
                         ry      = as.data.frame(!is.na(Xy_mis_ld))[, j], 
                         x       = as.matrix(Zm_ld[, -j]),
                         ls.meth = "ridge", ridge = ridge_trial)
  pdraw_hd <- .norm.draw(y       = Zm_hd[, j],
                         ry      = as.data.frame(!is.na(Xy_mis_hd))[, j], 
                         x       = as.matrix(Zm_hd[, -j]),
                         ls.meth = "ridge", ridge = ridge_trial)
  
  # Inside of .norm.draw w/ ls.meth = "ridge", this is what happens:
  # arguments
  x     = as.matrix(Zm_hd[, -j])
  ridge = ridge_trial
  
  # procedure
  xtx <- crossprod(x)
  pen <- ridge * diag(xtx)
  # if (length(pen) == 1) 
  #   pen <- matrix(pen)
  v <- solve(xtx + diag(pen))
  c <- t(y) %*% x %*% v
  r <- y - x %*% t(c)

  # Inspect posterior draws
  pdraw_ld$beta[1:49,]
  pdraw_ld$sigma
  
  pdraw_hd$beta[1:49,]
  pdraw_hd$sigma
  
# Interaction terms and number of predictors ------------------------------

  # Generaete some X predictors
  p <- 30 # number of variables
  X <- MASS::mvrnorm(1e2, rep(0, p), diag(p))
  
  # How many possible two-way interactions between these variables?
  n <- ncol(X)
  k <- 2 # two-way interactions 
  
  # Combination
  factorial(n) / (factorial(k)*factorial(n-k))
  
# Best Crossvalidation for Elastic net ------------------------------------

reps <- 25
lambda_val <- matrix(NA, ncol=3, nrow=reps)
alpha_val <- matrix(NA, ncol=3, nrow=reps)

for (i in 1:reps) {
  print(i)
  Xy <- genData(cond, parms)
  X <- Xy[, -201]
  y <- Xy$y
  
  # Version 1
  train_control <- trainControl(method = "cv",
                                number = 10,
                                # repeats = 5,
                                search = "random", # NEEDS CHECKING: what is this search method?
                                selectionFunction = "best",
                                verboseIter = FALSE)
  # CV
  el_cv <- train(y ~., 
                 data = cbind(y, X), 
                 method = "glmnet",
                 family = glmfam,
                 trControl = train_control,
                 tuneLength = 25
  )
  
  lambda_val[i, 1] <- as.numeric(el_cv$bestTune["lambda"])
  alpha_val[i, 1] <- as.numeric(el_cv$bestTune["alpha"])
  
  # Version 2
  train_control <- trainControl(method = "cv",
                                number = 10,
                                # repeats = 5,
                                selectionFunction = "best",
                                verboseIter = FALSE)
  # CV
  el_cv <- train(y ~., 
                 data = cbind(y, X), 
                 method = "glmnet",
                 family = glmfam,
                 trControl = train_control,
                 tuneLength = 25
  )
  lambda_val[i, 2] <- as.numeric(el_cv$bestTune["lambda"])
  alpha_val[i, 2] <- as.numeric(el_cv$bestTune["alpha"])
  
  # Version 3
  train_control <- trainControl(method = "cv",
                                number = 10,
                                repeats = 5,
                                selectionFunction = "best",
                                verboseIter = FALSE)
  # CV
  el_cv <- train(y ~., 
                 data = cbind(y, X), 
                 method = "glmnet",
                 family = glmfam,
                 trControl = train_control,
                 tuneLength = 25
  )
  lambda_val[i, 3] <- as.numeric(el_cv$bestTune["lambda"])
  alpha_val[i, 3] <- as.numeric(el_cv$bestTune["alpha"])
 
}

apply(lambda_val, 2, var)
apply(alpha_val, 2, var)


# Multicategory data gen check --------------------------------------------

# genMulticat check
X <- matrix(rnorm(1e5), ncol = 2) # some dataset
out <- genMulticat(K = 5, X)
y <- out$y
fit <- multinom(y ~ X)
list(est = round(coef(fit), 3),
     true = round(out$true_par,3))
# estimated and true parameters values are almost the same

# Correlation: single component structure --------------------------------------

# Define number of variables
p <- 50
rho <- seq(0, .9, by = .1)

# Storing objects
ngs <- NULL
cpve <- NULL

# Loop over conditions
for (r in seq_along(rho)) {
  # Define Sigma
  Sigma <- matrix(rho[r], p, p)
  diag(Sigma) <- 1

  # Sample data
  X <- MASS::mvrnorm(1e4, mu = rep(1, p), Sigma = Sigma)

  # PCA
  svd_X <- svd(X)

  # Compute cumulative proportion of explained variance
  cpve <- rbind(cpve, cumsum(prop.table(svd_X$d^2)))

  # Non-graphical solutions
  ngs <- rbind(ngs, nFactors::nScree(cor(X))$Components)
}

# Number of factors underlying data
cbind(
  rho = rho,
  ngs,
  npc.cpve = round(cpve[, 1:5], 1)
)

# Collinearity data check ------------------------------------------------------

rm(list = ls())
source("./init_general.R")
source("./exp1_init.R")

# Define experimental factor levels
p <- c(50) # c(50, 500) # number of variables
collinearity <- c(NA, seq(0.1, .9, by = .1))

# Create experimental conditions
conds <- expand.grid(
  p = p,
  collinearity = collinearity
)

# Increase sample size to make things more clear
parms$n <- 1e4

# Crate a place to store factor structures
storenScree <- list()
cpve <- NULL
pc1.loadings <- NULL
plots <- NULL
npcs_kpet <- NULL

# Define a plotting arrangement
par(mfrow = c(ceiling(sqrt(nrow(conds))), sqrt(nrow(conds))))

# Loop over the conditions
for (i in 1:nrow(conds)) {
  # Define the active condition
  cond <- conds[i, ]

  # Generate data
  Xy <- simData_exp1(cond, parms)

  # Scale data
  Xy <- scale(Xy)

  # Check correlation matrix
  print(round(cor(Xy)[c(1:13, 48:50), c(1:13, 48:50)], 1) * 100)

  # Select possible auxiliary data
  Xy_ma <- Xy[, -c(1:3, 6:8)]

  # 3-factor structure
  storenScree <- rbind(storenScree, nFactors::nScree(as.data.frame(Xy_ma))$Components)

  # PCA
  svd_X <- svd(Xy_ma)

  # Store the loadings
  pc1.loadings <- rbind(pc1.loadings, svd_X$v[, 4])

  # Compute CPVE for this run
  cpve_i <- cumsum(prop.table(svd_X$d^2))

  # Compute cumulative proportion of explained variance
  cpve <- rbind(cpve, cond = cpve_i)

  # Apply the decision rule used in the simulation study
  if (cpve_i[1] >= 0.5) {
    npcs_kpet <- c(npcs_kpet, 1)
  } else {
    npcs_kpet <- c(npcs_kpet, sum(cpve_i <= 0.5))
  }

  # Fit model regressing one variable on all the others
  model_all <- lm(z1 ~ ., data = as.data.frame(Xy))

  # Compute the VIFs
  vif_values <- car::vif(model_all)

  # Plot the vif values
  barplot(vif_values,
    main = paste0("VIF Values", " (cor = ", cond[, "collinearity"], ")"),
    horiz = TRUE,
    col = "steelblue"
  ) # create

  # And plot a reference line at 5 (VIF > 5 -> problem!)
  abline(v = 5, lwd = 3, lty = 2) # add vertical line at 5 as
  # If the value of VIF is :
  # - VIF < 1: no correlation
  # - 1 < VIF < 5, there is a moderate correlation
  # - VIF > 5: severe correlation
}

# Number of factors underlying data
cbind(
  collinearity = collinearity,
  storenScree,
  npcs50cpve = npcs_kpet,
  cpve = round(cpve[, 1:3], 2)*100
)

# Loadings
rownames(pc1.loadings) <- collinearity
colnames(pc1.loadings) <- colnames(Xy_ma)
round(abs(pc1.loadings), 1)[, c(1:4, 40:44)] * 10

# Test MI-PCA: run all collinearity values with 50% rule -----------------------

# > Prepare data ---------------------------------------------------------------

# Load results
out_MIPCA_colli <- readRDS(paste0("../output/", "exp1_2_simOut_20230421_1424.rds"))

# Review conditions run
out_MIPCA_colli$conds

# Model: means, variances, covariances (MLE estimates) per conditions
sem_res <- lapply(
  1:length(out_MIPCA_colli[[1]]),
  function(x) {
    res_sum(out_MIPCA_colli,
      model = "sem",
      condition = x
    )
  }
)

# Model: Linear regression per condition
lm_res <- lapply(
  1:length(out_MIPCA_colli[[1]]),
  function(x) {
    res_sum(out_MIPCA_colli,
      model = "lm",
      condition = x
    )
  }
)

# Shape for plot
output <- lapply(
  list(
    sem = sem_res,
    lm = lm_res
  ),
  function(x) {
    names(x) <- paste0("cond", seq_along(out_MIPCA_colli[[1]]))
    return(x)
  }
)
output$parms <- out_MIPCA_colli$parms
output$conds <- out_MIPCA_colli$conds

gg_out_sem <- plotwise(
  res = output,
  model = "sem",
  parPlot = list(
    Means = 1:6,
    Variances = 7:12,
    Covariances = 13:27
  ),
  item_group = c(1:3), # items in a group receiving miss values
  meth_compare = c(
    "MI_PCA",
    "MI_am",
    "CC",
    "GS"
  ),
  exp_factors = c("p", "collinearity")
)

# > Bias plot ------------------------------------------------------------------

# Bias or CIC?
x <- 1 # bias

# Which methods
methods_sel <- levels(gg_out_sem$methods)

# Which p condition
p_grep <- 500

# X breaks
xci_breaks <- sort(c(0, 10, 20, 50))

# Plot font
plot_text_family <- "sans"
plot_text_face <- "plain"
plot_text_size <- 9

# Make the plot
pf <- gg_out_sem %>%
    filter(
        analysis == unique(analysis)[x],
        grepl(p_grep, cond),
        variable %in% c("Min", "Mean", "Max"),
        methods %in% methods_sel
    ) %>%
    mutate(methods = fct_relabel(
        methods,
        str_replace,
        "-la", ""
    )) %>%
    # Main Plot
    ggplot(data = ., aes(
        y = methods,
        x = value,
        shape = variable
    )) +
    geom_point(size = 1.75) +
    geom_line(aes(group = methods),
        size = .25
    ) +
    # Grid
    facet_grid(
        rows = vars(factor(parm,
            levels = unique(parm)
        )),
        cols = vars(cond)
    ) +
    geom_vline(
        data = data.frame(
            xint = 10,
            analysis = "Percentage Relative Bias"
        ),
        linetype = "solid",
        size = .15,
        aes(xintercept = xint)
    ) +

    # Format
    scale_x_continuous(
        labels = xci_breaks,
        breaks = xci_breaks
    ) +
    scale_y_discrete(limits = rev) +
    scale_shape_manual(values = c("I", "I", "I")) +
    coord_cartesian(xlim = c(0, 50)) +
    labs( # title = label_parm[x],
        x = NULL,
        y = NULL,
        linetype = NULL,
        shape = NULL
    ) +
    theme(
        panel.background = element_rect(
            fill = NA,
            color = "gray"
        ),
        panel.grid.major = element_line(
            color = "gray",
            size = 0.15,
            linetype = 1
        ),
        legend.key = element_rect(
            colour = "gray",
            fill = NA,
            size = .15
        ),
        text = element_text(
            family = plot_text_family,
            face = plot_text_face,
            size = plot_text_size
        ),
        axis.ticks = element_blank(),
        legend.position = "none"
    )

pf

# > CIC plot -------------------------------------------------------------------

x <- 2 # Confidence intervals
methods_sel <- levels(gg_out_sem$methods)[1:2] # [-8]

# SE for threshold
ci_lvl <- .95
dt_reps <- 500
SEp <- sqrt(ci_lvl * (1 - ci_lvl) / dt_reps)
low_thr <- (.95 - SEp * 2) * 100
hig_thr <- (.95 + SEp * 2) * 100
vline_burton <- c(low_thr, hig_thr)
vline_vanBuu <- c(90, 99)
xci_breaks <- sort(c(80, vline_vanBuu, 95, round(vline_burton, 1), 100))

pf <- gg_out_sem %>%
  filter(
    analysis == unique(analysis)[x],
    grepl(p_grep, cond),
    variable %in% c("Min", "Mean", "Max"),
    methods %in% methods_sel
  ) %>%
  # Drop pm = ** as we are plotting only one value for this condition
  mutate(cond = fct_relabel(
    cond,
    str_replace,
    " pm = [0-9]\\.[0-9]", ""
  )) %>%
  mutate(methods = fct_relabel(
    methods,
    str_replace,
    "-la", ""
  )) %>%
  # Main Plot
  ggplot(data = ., aes(
    y = methods,
    x = value,
    shape = variable
  )) +
  geom_point(size = 1.75, show.legend = FALSE) +
  geom_line(aes(group = methods),
    size = .25
  ) +
  # Grid
  facet_grid(
    rows = vars(factor(parm,
      levels = unique(parm)
    )),
    cols = vars(cond)
  ) +
  geom_vline(
    data = data.frame(
      xint = 95,
      analysis = "CI coverage"
    ),
    linetype = "solid",
    size = .15,
    aes(
      xintercept = xint,
      lty = paste0("nominal level")
    )
  ) +

  # Format
  scale_y_discrete(limits = rev) +
  scale_x_continuous(
    labels = as.character(round(xci_breaks / 100, 2)),
    breaks = xci_breaks
  ) +
  coord_cartesian(xlim = c(min(xci_breaks), max(xci_breaks))) +
  scale_shape_manual(values = c("I", "I", "I")) +
  labs( # title = label_parm[x],
    x = NULL,
    y = NULL,
    linetype = NULL,
    shape = NULL
  ) +
  theme(
    panel.background = element_rect(
      fill = NA,
      color = "gray"
    ),
    panel.grid.major = element_line(
      color = "gray",
      size = 0.175,
      linetype = 1
    ),
    axis.ticks = element_blank(),
    legend.key = element_rect(
      colour = "gray",
      fill = NA,
      size = .15
    ),
    legend.position = "bottom",
    text = element_text(
      family = plot_text_family,
      face = plot_text_face,
      size = plot_text_size
    ),
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    )
  )

pf
