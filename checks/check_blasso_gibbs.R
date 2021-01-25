### Title:    Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-06-08
### Notes:    In your code you are using the blasso.vs funcion to obtain 
###           one sample at the time from the coditional posteriors within
###           the larger gibbs sampler for multiple imputation. This scripts
###           makes sure that running the gibbs sampler in this way is 
###           equivalent to drawing from each conditional distirbtuion.

rm(list=ls())
source("./init_general.R")
set.seed(20200617)

# Get data
N <- 500
p <- 10

# X
X <- rmvnorm(N, 
             mean = rep(0, p),
             sigma = AR1(p, .2))
X_design <- cbind(1, X)

# True parms
B <- c(b0 = 0,
       c(1.5, 3), 
       rep(.001, p-2))

# Y
y <- X_design %*% B + rnorm(N)

# Gibbs sampler details ---------------------------------------------------

# Iterations
iters <- 1e4
bin <- 1e4

# Initial values
B0 <- rep(0, ncol(X))
tau0 <- 1
sigma20 <- 1
phi0 <- .5

# Straightforward use of blasso -------------------------------------------

dir_blasso <- blasso::blasso.vs(y, X, 
                                iters = iters,
                                burn  = bin, 
                                beta  = B0, 
                                beta.prior = "scaled", 
                                sig2 = sigma20, sig2prior = c(a = .1 , b = .1),
                                tau = tau0, tauprior=c(r = 1, s = 1),
                                phi = phi0, phiprior=c(g = 1, h = 1),
                                fixsig = FALSE, fixtau = FALSE, fixphi = FALSE)

round(colMeans(dir_blasso$beta), 3)

plot(1:iters, dir_blasso$beta[, 1], type = "l",
     main = paste0("b", 1:p)[1],
     xlab = "iterations",
     ylab = "")
plot(density(dir_blasso$beta[, 1]),
     main = paste0("b", 1:p)[1],
     xlab = "coef values",
     ylab = "density")

par(mfrow=c(3,2))
lapply(1:ncol(dir_blasso$beta), 
       function(x) {
         plot(1:iters, dir_blasso$beta[, x], type = "l",
              main = paste0("b", 1:p)[x],
              xlab = "iterations",
              ylab = "")
         plot(density(dir_blasso$beta[, x]),
              main = paste0("b", 1:p)[x],
              xlab = "coef values",
              ylab = "density")
       }
)

# Single draw repeated ----------------------------------------------------

# regression coefficients
Beta.out <- matrix(NA, nrow = bin+iters, ncol = ncol(X))
  Beta.out[1, ] <- 0

# Error variance
Sig.out <- matrix(NA, nrow = bin+iters, ncol = 1)
  Sig.out[1] <- 1
  
# tau parameter in Hans 2009 and 10 blasso model
Tau.out <- matrix(NA, nrow = bin+iters, ncol = 1)
  Tau.out[1] <- 1

# phi parameter in Hans 2010 blasso model
Phi.out <- matrix(NA, nrow = bin+iters, ncol = 1)
  Phi.out[1] <- .5

for (m in 2:(bin+iters)) {
  beta_m  <- Beta.out[m-1,]
  sigma_m <- sqrt(Sig.out[m-1])
  tam_m   <- Tau.out[m-1]
  phi_m   <- Phi.out[m-1]

  pdraw <- blasso::blasso.vs(Y = y, X = X,
                             iters = 1,
                             burn  = 0,
                             beta  = beta_m, 
                             beta.prior = "scaled", 
                             sig2  = sigma_m, sig2prior = c(a = .1, b = .1),
                             tau   = tam_m,   tauprior  = c(r = .01, s = .01), 
                             phi   = phi_m,   phiprior  = c(h = 1, g = 1),
                             fixsig = FALSE, fixtau = FALSE, fixphi = FALSE,
                             noisy = FALSE)
    
  # Store results
  Beta.out[m, ] <- as.vector(pdraw$beta)
  Sig.out[m]    <- pdraw$sig2
  Tau.out[m]    <- pdraw$tau
  Phi.out[m]    <- pdraw$phi
}
  
  post_beta <- Beta.out[-c(1:bin), ]
  post_sig2 <- Sig.out[-c(1:bin), ]
  post_tau <- Tau.out[-c(1:bin), ]
  post_phi <- Phi.out[-c(1:bin), ]
  
  par(mfrow=c(3,2))
  lapply(1:ncol(X), 
         function(x) {
           plot(1:iters, post_beta[, x], type = "l",
                main = paste0("b", 1:p)[x],
                xlab = "iterations",
                ylab = "")
           plot(density(post_beta[, x]),
                main = paste0("b", 1:p)[x],
                xlab = "coef values",
                ylab = "density")
         }
  )
  

# Compare posteriors 1 on 1 -----------------------------------------------

# ---------------------- #
# Regression coefficients
# ---------------------- #
  
  # Posterior summaries
  data.frame(direct = round(colMeans(dir_blasso$beta), 3),
             man_it = round(colMeans(post_beta), 3))
  
  # Traceplot and densities
  x <- 1 # select 1 beta
  par(mfrow=c(1,3))
  
  # Direct use
  plot(1:iters, dir_blasso$beta[, x], type = "l",
       main = paste0("b", 1:p)[x],
       xlab = "iterations",
       ylab = "")
  
  # Use within iterations
  plot(1:iters, post_beta[, x], type = "l",
       main = paste0("b", 1:p)[x],
       xlab = "iterations",
       ylab = "")
  
  # Direct use vs w/in iterations
  plot(density(dir_blasso$beta[, x]),
       main = paste0("b", 1:p)[x],
       xlab = "coef values",
       ylab = "density", lwd = 2)
  lines(density(post_beta[, x]),
        col = "blue", lty = 2, lwd = 2)
  
# ---------------------- #
# Sigma
# ---------------------- #
  
  plot(density(dir_blasso$sig2),
       main = "sigma2",
       xlab = "coef values",
       ylab = "density", lwd = 2)
  lines(density(post_sig2),
        col = "blue", lty = 2, lwd = 2)
  
# ---------------------- #
# Tau
# ---------------------- #
  
  plot(density(dir_blasso$tau),
       main = "tau",
       xlab = "coef values",
       ylab = "density", lwd = 2)
  lines(density(post_tau),
        col = "blue", lty = 2, lwd = 2)
  
# ---------------------- #
# Phi
# ---------------------- #
  
  plot(density(dir_blasso$phi),
       main = "phi",
       xlab = "coef values",
       ylab = "density", lwd = 2)
  lines(density(post_phi),
        col = "blue", lty = 2, lwd = 2)