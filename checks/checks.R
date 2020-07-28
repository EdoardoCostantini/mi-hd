### Title:    Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-05-19

rm(list=ls())
source("./init_general.R")
source("./init_exp1.R")
source("./init_exp3.R")

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
source("./init_exp1.R")

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

# Paramters of interest ---------------------------------------------------

rm(list=ls())
source("./init.R")
set.seed(1234)

cond <- conds[3,]

reps <- 1e2
full_store <- matrix(NA, nrow = reps, ncol = length(parms$z_m_id))
miss_store <- matrix(NA, nrow = reps, ncol = length(parms$z_m_id))

full_cov_sum <- matrix(0, nrow = length(parms$z_m_id), ncol = length(parms$z_m_id))
miss_cov_sum <- matrix(0, nrow = length(parms$z_m_id), ncol = length(parms$z_m_id))

full_cor_sum <- matrix(0, nrow = length(parms$z_m_id), ncol = length(parms$z_m_id))
miss_cor_sum <- matrix(0, nrow = length(parms$z_m_id), ncol = length(parms$z_m_id))

for (r in 1:reps) {
  print(r/reps*100)
  dt_full <- genDt_mvn(cond, parms)
  dt_mis  <- imposeMiss(dt_full, parms, cond)
  full_store[r, ] <- colMeans(dt_full)[parms$z_m_id]
  miss_store[r, ] <- colMeans(dt_mis, na.rm = TRUE)[parms$z_m_id]
  
  full_cov_sum <- full_cov_sum + cov(dt_full[, parms$z_m_id])
  miss_cov_sum <- miss_cov_sum + cov(dt_mis[, parms$z_m_id], use = "pairwise.complete.obs")
  
  full_cor_sum <- full_cor_sum + cor(dt_full[, parms$z_m_id])
  miss_cor_sum <- miss_cor_sum + cor(dt_mis[, parms$z_m_id], use = "pairwise.complete.obs")
}

# Effect of missingness on analysis

# Means

  round(
    cbind(colMeans(full_store),
          colMeans(miss_store)),
    3)
  
  full_mv <- colMeans(full_store)
  miss_mv <- colMeans(miss_store)
  
  # Bias in terms of percentage of true value
  round((full_mv - miss_mv)/full_mv*100, 0)

# COVARIANCE
  
  full_cov <- full_cov_sum/reps
  miss_cov <- miss_cov_sum/reps
  lapply(list(full_cov, miss_cov), round, 2)
  # Bias percent points (percentage of true value)
  bias_pg <- round((full_cov - miss_cov)/full_cov*100, 0)
  bias_pg
  pm_50 <- bias_pg[lower.tri(bias_pg)]
  data.frame(pm_50 = pm_50, pm_500 = pm_500)
  
  mean(diag(bias_pg)) # mean variances bias %
  mean(bias_pg[lower.tri(bias_pg)]) # mean covariances bias %

# CORRELATION
  
  full_cor <- full_cor_sum/reps
  miss_cor <- miss_cor_sum/reps
  lapply(list(full_cor, miss_cor), round, 2)
  # Bias in terms of percentage of true value
  round((full_cor - miss_cor)/full_cor*100, 0)
  mean(unique(round((full_cor - miss_cor)/full_cor*100, 0)))
  

# Ridge Penalty choice ----------------------------------------------------

  rm(list=ls())
  source("./init_general.R")
  source("./init_exp1.R")
  
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
  

# Data Generation Latent Structure ----------------------------------------

  rm(list=ls())
  source("./init_general.R")
  source("./init_exp3.R")
  
## Works for all conditions ##
  
  all_conds_data <- lapply(1:nrow(conds), function(x){
    X <- simData_exp3(parms, conds[x,])
    Xy_mis <- imposeMiss_lv(X, parms, conds[x, ])
    return(Xy_mis)
  })
  length(all_conds_data)
  
## On average, items have mean 0, variance 1 ##
  
  set.seed(20200727)
  store <- matrix(0, nrow = parms$n_it*conds[1,]$lv, ncol = 2)
  for (i in 1:5e3) {
    X <- simData_exp3(parms, conds[1,])
    store <- store + data.frame(mu = colMeans(X$dat), 
                                var = apply(X$dat, 2, var))
  }
  round(store/5e3, 3)
  
## check estimates on raw original data ##
  
  # Deifne some CFA mode to test
  CFA_model <- '
    # Measurement Model
    lv1 =~ z1 + z2 + z3 + z4 + z5
    lv2 =~ z6 + z7 + z8 + z9 + z10
  '
  
  set.seed(20200727)
  X <- simData_exp3(parms, conds[1,])
  
  # Scaled data
  Xs <- X$dat * sqrt(5) + 5

  Xs <- sapply(X$dat, function(x){
    x * sqrt(5)/sd(x) + 5
  })
  colMeans(X$dat)
  apply(X$dat, 2, var)
  colMeans(Xs)
  apply(Xs, 2, var)
  
  # Fit models
  fit   <- cfa(CFA_model, data = X$dat, std.lv = TRUE)
  fit_v <- cfa(CFA_model, data = X$dat)
  fit_s <- cfa(CFA_model, data = Xs, std.lv = TRUE)
  fit_s <- cfa(CFA_model, data = scale(Xs), std.lv = TRUE)
  
  # Measured
  round(X$Lambda, 3)  # factor loadings
  diag(X$Theta)[1:10] # items error variances
  
  # Setting latent variables variances to 1 standardizes factor scores already
  parameterEstimates(fit, 
                     se = F, zstat = F, pvalue = F, ci = F,
                     standardized = TRUE)
  parameterEstimates(fit_v, 
                     se = F, zstat = F, pvalue = F, ci = F,
                     standardized = TRUE)
  parameterEstimates(fit_s,
                     se = F, zstat = F, pvalue = F, ci = F,
                     standardized = TRUE)
  
  # Compare
  round(X$Lambda, 3)[1,1] / round(X$Lambda, 3)[2,1]
  parameterEstimates(fit)$est[1] /parameterEstimates(fit)$est[2]
  parameterEstimates(fit_s, standardized = TRUE)$est[1] / parameterEstimates(fit_s, standardized = TRUE)$est[2]

## Low / high factor loadings
  set.seed(20200727)
  X_high <- simData_exp3(parms, conds[1,])
  set.seed(20200727)
  X_low  <- simData_exp3(parms, conds[5,])
  
  fits <- lapply(list(X_high$dat, X_low$dat), 
                 cfa, model = CFA_model, std.lv = TRUE)
  lapply(fits, parameterEstimates, standardized = TRUE,
         se = F, zstat = F, pvalue = F, ci = F)
  round(X_high$Lambda, 3)
  round(X_low$Lambda, 3)

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