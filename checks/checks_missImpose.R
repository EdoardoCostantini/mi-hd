### Title:    Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-07-08

# Missing data pm ---------------------------------------------------------
# Does the proportion of missing data imposed work?

rm(list=ls())
source("./init_general.R")
source("./init_exp1.R")

reps <- 1e3
cond <- conds[1,]
pm_store <- matrix(NA, nrow = reps, ncol = cond$p)

for (r in 1:reps) {
  dt_full <- simData_exp1(cond, parms)
  dt_mis  <- imposeMiss(dt_full, parms, cond)
  pm_store[r, ] <- colMeans(is.na(dt_mis))
}

round(colMeans(pm_store), 3)
round(colMeans(pm_store), 1) == cond$pm
# All variables with missing values have on average
# the correct proportion of missing cases.

# Missing Mechanism Understanding -----------------------------------------
# How deos the simMissngness.R script work?

# > Set up ####

  rm(list=ls())
  source("./init_general.R")
  source("./init_exp1.R")
  
  ## gen data according to condition ##
  cond <- conds[1, ]
  Xy <- simData_exp1(cond, parms)
  dim(Xy)
  data  = Xy
  preds = parms$rm_x[1, ]

# > Internals ####
  # Internals
  # Type of missingnes (where in the target variable distribution)
    type = "center"
    
  # Eta: Standardized linear predictors
    set.seed(1234)
    preds = c(2)
    beta = rep(3, length(preds)) 
      # as long as they are the same value, eta will be exactly the same
      # these values define the relative importance of each predictor
    data = Xy
    eta = scale(
      as.matrix(data[ , preds]) %*% matrix(beta)
    )
    
    # Check 
    c(mean(as.matrix(data[ , preds]) %*% matrix(beta)), 
      sd(as.matrix(data[ , preds]) %*% matrix(beta)))
    c(mean(eta), sd(eta))
    plot(density(eta)) 
      # the best check: stays exactly the same as long as beta values keep
      # same ratio (same value for all, means they all weight the same,
      # first one is double in size than the others, means it weights 
      # twice as much as each of the others)
    
  # Proportion of missing cases (on target variable)
    pm = .3
  
  # Optimize the offset:
    offset <- optOffset(pm = pm, eta = eta, type = type)$minimum
  
  # > Step 1: "fOffset" ####  
  # Create the objective function for logistic regression offset
    
    fOffset(offset, eta, pm, type)
    # Creates some offsetting value that will be used later
    # This is just hte function that given certain parameters returns a value of the
    # ofset, not the one that actually retunr the offset values to be used.
    fOffset <- function(offset, eta, pm, type) {
      f <- switch(type,
                  center = offset - abs(eta), # puts missings in the center
                  tails  = abs(eta) - offset, # puts missings in the tails
                  offset + eta                # puts missings in one tail (high or low, 
                                              # sign not important at this stage)
      )
      
      (mean(plogis(f)) - pm)^2
    }

# > Step 2: "optOffset" ####  
  
  ## Optimize the logistic regression offset for a given value of the linear
  ## predictor (eta) to get a desired percent missing (pm).

  # Max iterations for optimazier
  maxIter = 10
  
  # Tolarance for convergence of optimazier
  tol = c(0.1, 0.001)
  
  # obtain actual offset value
  offset <- optOffset(pm = pm, eta = eta, type = type)$minimum

  # How does it work? Basically just uses omptimize to find 
  # the value of offset that minimizes the objective function
  # discussed above.
  
  # optOffset <- function(pm, eta, type, tol = c(0.1, 0.001), maxIter = 10) {
  #   for(k in 1 : maxIter) {
      k <- 1
  #     ## Define the search range:
      int <- k * range(eta)
      
      ## Optimize the objective over 'int':
      out <- optimize(f        = fOffset,
                      interval = int,
                      eta      = eta,
                      pm       = pm,
                      type     = type)
      
  #     ## Are we far enough from the boundary?
  #     dist   <- out$minimum - int
  #     check1 <- all(abs(dist) > tol[1] * diff(int))
  #     
  #     ## Are we within tolerance?
  #     check2 <- out$objective < tol[2]
  #     
  #     if(check1 & check2) break
  #   }
  #   ## Did we fail?
  #   if(!check1 | ! check2) stop("I could not optimize this function.")
  #   
  #   out
  # }
    
# > Step 3: "simMissingness" ####
    
  # Actually generates the response vector based on the standard 
  # logistic CDF.
    
  ## Simulate a nonresponse vector:
  # simMissingness <- function(pm,
  #                            data,
  #                            preds = colnames(data),
  #                            type  = "high",
  #                            beta  = NULL)
  # {
  #   ## Define a trivial slope vector, if necessary:
  #   if(is.null(beta)) beta <- rep(1.0, length(preds))
  #   
  #   ## Compute (and standardize) the linear predictor:
  #   eta <- scale(
  #     as.matrix(data[ , preds]) %*% matrix(beta)
  #   )
  #   
  #   ## Optimize the offset:
  #   offset <- optOffset(pm = pm, eta = eta, type = type)$minimum
    
    ## Compute the probabilities of nonresponse:
    # Based on the object "type" decide what to do with the
    # offset object that was created with the optOffset function
    probs <- plogis(
      switch(type,
             high   = offset + eta,
             low    = offset - eta,
             center = offset - abs(eta),
             tails  = abs(eta) - offset
      )
    )

    ## Return a logical nonresponse vector:
    as.logical(rbinom(n = length(eta), size = 1, prob = probs))
    
  # }

# > Visualize ####
    set.seed(20200708)
    x <- sort(rlogis(1e5,
                     location = 0, # standard logistic
                     scale = 1))
    id <- order(eta)[5]
      # some observation in my dataset that is in the low part of the data
    id <- order(eta)[195]
      # some observation in my dataset that is in the high part of the data
    
  # Density of fully observed vs var with missings
    
  # Probabilities of missings for a given id, w/ and w/out offset
    sort(
      round(
        c(raw    = plogis(eta[id]),
          high   = plogis(offset + eta[id]), 
          low    = plogis(offset - eta[id]), 
          center = plogis(offset - abs(eta[id])),
          tails  = plogis(abs(eta[id]) - offset)),
        3
      )
    )
    
  # Plot PDF and CDF
    par(mfrow = c(1, 2))
    plot(density(x), 
         main = "PDF", xlab = "x")
      points(eta[id], dlogis(eta[id]))
      abline(v = eta[id], lty = 2, col = "gray")
      text(eta[id]+2.5, dlogis(eta[id]), labels = paste0("id = ", id), col = "gray")
    plot(x, plogis(x), 
         main = "CDF", ylab = "Probability", type = "l")
    
  # Compute Probabilities of nonresponse w/out offset
    probs <- plogis(eta[id])
      points(eta[id], probs)
      abline(h = probs, lty = 2, col = "gray")
      text(-7.5, probs+.025, 
           labels = paste0("p(NA) = ", round(probs, 2)), 
           col = "gray", cex = .75)
      
    
  # Compute Probabilities of nonresponse w/ offset low tail
    probs <- plogis(offset - eta[id])
      points(eta[id], probs)
      abline(h = probs, lty = 2, col = "gray")
      text(-7.5, probs+.025, 
           labels = paste0("p(NA) = ", round(probs, 2)), 
           col = "red", cex = .75)
      
  # Compute Probabilities of nonresponse w/ offset high tail
    probs <- plogis(offset + eta[id])
      points(eta[id], probs)
      abline(h = probs, lty = 2, col = "gray")
      text(-7.5, probs+.025, 
           labels = paste0("p(NA) = ", round(probs, 2)), 
           col = "blue", cex = .75)
      
  # Compute Probabilities of nonresponse w/ offset high tail
    probs <- plogis(offset - abs(eta[id]))
      points(eta[id], probs)
      abline(h = probs, lty = 2, col = "gray")
      text(-7.5, probs+.025, 
           labels = paste0("p(NA) = ", round(probs, 2)), 
           col = "blue", cex = .75)
      
  # Compute Probabilities of nonresponse w/ offset high tail
    probs <- plogis(abs(eta[id]) - offset)
      points(eta[id], probs)
      abline(h = probs, lty = 2, col = "gray")
      text(-7.5, probs+.025, 
           labels = paste0("p(NA) = ", round(probs, 2)), 
           col = "red", cex = .75)
      

# Density Comparison ------------------------------------------------------
# Compare fully observed variable and different versions of one with 
# missing values imposed with offsetting. In this run, you can see
# how for different type of imposition (low, high, center, tails)
# the density of the observed part of variable is different from the
# fully observed one.
      
  rm(list=ls())
  source("./init_general.R")
  source("./init_exp1.R")
      
  parms$n <- 5e4 # increased to make patterns clearer
  
  # Gen data
  set.seed(20200805)
  cond <- conds[1, ] # actual pm is defined below
  Xy <- simData_exp1(cond, parms)
  Xy_mis <- Xy
  
  # Plot density of fully observed variable (black solid)
  plot(density(Xy[, 1]),
       main = "Effects of offsetting missImpose",
       ylim = c(0, max(density(Xy[, 1])$y) + sd(density(Xy[, 1])$y)),
       lwd = 3
  )
  
  # Chose which methods to plot
  range <- 1:4
  types <- c("low", "high", "center", "tails")[range]
  color <- c("red", "blue", "darkorange", "darkgray")[range]
  
  # Plot density of observed part under the different methods
  lapply(range, function(x){
    for (i in 1:50) {
      nR <- simMissingness(data  = Xy, 
                           preds = c(2,3,4,5), # more predictors accentuate effect
                           #beta  = 1, # any number gives same result
                                       # if all predictors have same importance
                           pm    = .8, # larger makes pattern clearer
                           type  = types[x])
      lines(density(Xy_mis[!nR, 1], na.rm = TRUE), 
            col = color[x], lty = 2, lwd = .5)
    }
  })
  legend(-2, .25,
         legend = types,
         col    = color, 
         lwd = c(2, 2, 2, 2),
         lty = c(2, 2, 2, 2), 
         cex = 1)
      
# Effect of missingness on analysis ---------------------------------------

# > Set up ####

  rm(list=ls())
  source("./init_general.R")
  source("./init_exp1.R")
      
  parms$auxWts <- rep(1, 4)
  parms$auxWts <- rep(4, 4) 
    # as long as they are equal among themselves,
    # it does not matter which number you are using.

# Loop data generation
  set.seed(1234)
  reps <- 5e2 # how many data repplications
  cond <- conds[3, ] # chosen conditions

  # Storing objects
  full_store <- matrix(NA, nrow = reps, ncol = length(parms$z_m_id))
  miss_store <- matrix(NA, nrow = reps, ncol = length(parms$z_m_id))
  
  full_cov_sum <- matrix(0, nrow = length(parms$z_m_id), ncol = length(parms$z_m_id))
  miss_cov_sum <- matrix(0, nrow = length(parms$z_m_id), ncol = length(parms$z_m_id))
  
  full_cor_sum <- matrix(0, nrow = length(parms$z_m_id), ncol = length(parms$z_m_id))
  miss_cor_sum <- matrix(0, nrow = length(parms$z_m_id), ncol = length(parms$z_m_id))

for (r in 1:reps) {
  print(r/reps*100)
  dt_full <- simData_exp1(cond, parms)
  dt_mis  <- imposeMiss(dt_full, parms, cond)
  full_store[r, ] <- colMeans(dt_full)[parms$z_m_id]
  miss_store[r, ] <- colMeans(dt_mis, na.rm = TRUE)[parms$z_m_id]
  
  full_cov_sum <- full_cov_sum + cov(dt_full[, parms$z_m_id])
  miss_cov_sum <- miss_cov_sum + cov(dt_mis[, parms$z_m_id], 
                                     use = "complete.obs")
  
  full_cor_sum <- full_cor_sum + cor(dt_full[, parms$z_m_id])
  miss_cor_sum <- miss_cor_sum + cor(dt_mis[, parms$z_m_id], 
                                     use = "complete.obs")
}

# Means
  means <- as.data.frame(sapply(list(GS = full_store,
                                     CC = miss_store), colMeans))
  # Bias in percent of reference value
  mu_bias_p <- round((means$GS - means$CC)/means$GS*100, 0)

# Variances and covariances
  vcovs <- lapply(list(GS = full_cov_sum/reps, 
                       CC = miss_cov_sum/reps), round, 2)
  vcovs
  
# Bias percent points (percentage of true value)
  bias_pg <- round((vcovs$GS - vcovs$CC)/vcovs$GS*100, 0)
  # All
  bias_pg
  
# Summary
  # Means
  mu_bias_p
  # Variances
  diag(bias_pg)
  mean(diag(bias_pg)) # mean variances bias %
  # Covariances
  bias_pg[lower.tri(bias_pg)]
  mean(bias_pg[lower.tri(bias_pg)]) # mean covariances bias %

# Correlation
  cors <- lapply(list(GS = full_cor_sum/reps, 
                       CC = miss_cor_sum/reps), round, 2)
  cors
  # Bias in terms of percentage of true value
  cor_bias <- round((cors$GS - cors$CC)/cors$GS*100, 0)
  # All
  cor_bias

      