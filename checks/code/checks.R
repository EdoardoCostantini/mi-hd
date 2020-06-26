### Title:    Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-05-19

source("./init.R")
source("./functions.R")

# Missing data pm ---------------------------------------------------------
# Does the proportion of missing data imposed work?

rm(list=ls())
source("./init.R")

reps <- 1e3
cond <- conds[1,]
pm_store <- matrix(NA, nrow = reps, ncol = cond$p)

for (r in 1:reps) {
  dt_full <- genDt_mvn(cond, parms)
  dt_mis  <- imposeMiss(dt_full, parms, cond)
  pm_store[r, ] <- colMeans(is.na(dt_mis))
}

round(colMeans(pm_store), 1) == cond$pm
# All variables with missing values have on average
# the correct proportion of missing cases.


# -------------------------------------------------------------------------

# SEM Model non convergence when few iterations
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


# Paramters of interest ---------------------------------------------------

rm(list=ls())
source("./init.R")
set.seed(1234)

cond <- conds[1,]

reps <- 1e3
full_store <- matrix(NA, nrow = reps, ncol = length(parms$z_m_id))
miss_store <- matrix(NA, nrow = reps, ncol = length(parms$z_m_id))

full_cov_sum <- matrix(0, nrow = length(parms$z_m_id), ncol = length(parms$z_m_id))
miss_cov_sum <- matrix(0, nrow = length(parms$z_m_id), ncol = length(parms$z_m_id))

full_cor_sum <- matrix(0, nrow = length(parms$z_m_id), ncol = length(parms$z_m_id))
miss_cor_sum <- matrix(0, nrow = length(parms$z_m_id), ncol = length(parms$z_m_id))

for (r in 1:reps) {
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
round(
  cbind(colMeans(full_store),
        colMeans(miss_store)),
  3)

full_mv <- colMeans(full_store)
miss_mv <- colMeans(miss_store)

# Bias in terms of percentage of true value
  round((full_mv - miss_mv)/full_mv*100, 0)

full_cov <- full_cov_sum/reps
miss_cov <- miss_cov_sum/reps
lapply(list(full_cov, miss_cov), round, 2)

# Bias percent points (percentage of true value)
  round((full_cov - miss_cov)/full_cov*100, 0)
  mean(unique(round((full_cov - miss_cov)/full_cov*100, 0)))

full_cor <- full_cor_sum/reps
miss_cor <- miss_cor_sum/reps
lapply(list(full_cor, miss_cor), round, 2)
# Bias in terms of percentage of true value
  round((full_cor - miss_cor)/full_cor*100, 0)
  mean(unique(round((full_cor - miss_cor)/full_cor*100, 0)))

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