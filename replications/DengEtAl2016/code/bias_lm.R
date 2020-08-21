
# Notes: CLosest attempt to replicate DengEtAl and ZaoLong GS and CC reuslts
#        When pp = 1, replicaiton of ZaoLong 2016 is successful within 
#        rounding error. When pp = 2, replication of Deng Et Al is unsuccessful
#        with bias that is within rounding error but se and sd very different
#        from reported in paper. Still I do get bias, whihc is what I need.

library("mvtnorm")
library("CVTuningCov")

rm(list=ls())

# Paper
pp <- 2 # 1 = ZaoLong, 2 = DengEtAl

nobs <- 1e2
p <- 200
reps <- 5e2

# Fully observed Xs
rho   <- c(0.5, .1)[pp]

# Var w/ miss model
n_mis <- c(1, 3)[pp] # number of variables with missings
s2_zm <- c(1, 16)[pp] # error vairance for var w/ miss
Z_s  <- list(c(2, 3, 50, 51),
             c(4, 5, 50, 51))[[pp]]
a <- rep(1, length(Z_s))

# Y model
s2_y <- c(3, 36)[pp] # error vairance for y var
y_mod_x <- list(1:3,
                1:5)[[pp]]
  b <- rep(1, length(y_mod_x)+1)
  
# Response model
rm_x <- c(2, 3) # which columns of x are predictors of miss
rm_x <- list(c(2, 3),
             data.frame(4:5,
                        c(4, 51),
                        51:52)
             )[[pp]]# which columns of x are predictors of miss
rm_coefs <- c(-1, -.1, 2, -2) # weights

# Storing objects
GS_b <- matrix(NA, nrow = reps, ncol = length(y_mod_x))
CC_b <- matrix(NA, nrow = reps, ncol = length(y_mod_x))
GS_se <- matrix(NA, nrow = reps, ncol = length(y_mod_x))
CC_se <- matrix(NA, nrow = reps, ncol = length(y_mod_x))
pm <- matrix(NA, nrow = reps, ncol = n_mis)

# set.seed(1234)

for (r in 1:reps) {
  Zf <- mvtnorm::rmvnorm(n     = nobs, 
                         mean  = rep(0, p-n_mis), 
                         sigma = AR1(p-n_mis, rho = rho) )
  Zm <- sapply(1:n_mis, 
               function(s){
                 rnorm(n = nobs, 
                       mean  = 1+Zf[, Z_s] %*% a, 
                       sd = sqrt(s2_zm) )
               })
  Zf[, 1:n_mis] <- Zm
  eps <- rnorm(nobs, 0, sd = sqrt(s2_y))
  y <- cbind(1, Zf[, y_mod_x]) %*% b + eps
  
  # Miss impose
  Xm <- Zf
  for (i in 1:n_mis) {
    if(pp == 1){
      logit_miss <- cbind(1, Xm[, rm_x], y) %*% rm_coefs
    } else {
      logit_miss <- cbind(1, Xm[, rm_x[[i]]], y) %*% rm_coefs
    }
    prb_miss   <- exp(logit_miss)/(1+exp(logit_miss))
    nR         <- rbinom(nobs, 1, prb_miss) == 1
    Xm[nR, i]  <- NA
  }
  
  datmis <- as.data.frame(cbind(y, Xm))
  O <- sapply(datmis, is.na)
  if(n_mis > 1){
    pm[r,]<- colMeans(O[, (1:n_mis)+1])
  } else {
    pm[r,]<- mean(O[, (1:n_mis)+1])
  }
  ind_cc <- rowSums(O) == 0
    
  lm_GS <- lm(y ~ Zf[, y_mod_x])
  lm_CC <- lm(y[ind_cc] ~ Xm[ind_cc, y_mod_x])
  
  GS_b[r,] <- coef(lm_GS)[-1]
  GS_se[r,]<- summary(lm_GS)$coefficients[-1, 2]
  CC_b[r,] <- coef(lm_CC)[-1]
  CC_se[r,]<- summary(lm_CC)$coefficients[-1, 2]
}

colMeans(pm)

bias = data.frame(GS=abs(colMeans(GS_b) - 1),
                  CC=abs(colMeans(CC_b) - 1))

se = data.frame(GS=abs(colMeans(GS_se)),
                CC=abs(colMeans(CC_se)))

sd = data.frame(GS=apply(GS_b, 2, sd),
                CC=apply(CC_b, 2, sd))

data.frame(t(round(bias[1, ], 3)),
           t(round(se[1, ],3)),
           t(round(sd[1, ], 3)))

round(
  data.frame(bias = as.vector(t(bias[1,])),
             se = as.vector(t(se[1,])),
             sd = as.vector(t(sd[1,])),
             row.names = c("GS", "CC")),
  3
)


t(round(bias[1:3, ], 3))
t(round(se[1:3, ], 3))
t(round(sd[1:3, ], 3))

