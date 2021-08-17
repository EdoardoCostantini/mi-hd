### Title:    Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-05-19

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
# estiamted and true parameters values are almost the same