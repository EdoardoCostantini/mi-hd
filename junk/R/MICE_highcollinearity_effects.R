# MICE imputation with collinearity issues
  # Simulation to investigate convergence problems in MICE algorithm
  # when strong correlation is present in multivariate dataset
  # This was performed on the lines of the simulation performed by
  # Buuren 2018 (p. 129, see code file under header "slowplot, echo=FALSE, fig.height=5")

library(mice)

# Functions ---------------------------------------------------------------

# Create a list of different correlation matrix
create_cor <- function(corr_coef){
  cormatrix <- NULL
  for (i in 1:length(corr_coef)) {
    cormatrix[[i]] <- matrix(c(1.0, 0.9, 0.9,
                             0.9, 1.0, corr_coef[i],
                             0.9, corr_coef[i], 1.0), nrow = 3)
  }
  return(cormatrix)
}
# Create dataset (fixed missing data pattern and proportion of missings)
generate <- function(n = c(1000, 4500, 4500, 0),
                     cor = matrix(c(1.0, 0.9, 0.9,
                                    0.9, 1.0, 0.7,
                                    0.9, 0.7, 1.0), nrow = 3)) {
  require(MASS)
  nt <- sum(n)
  cs <- cumsum(n)
  data <- mvrnorm(nt, mu = rep(0,3), Sigma = cor)
  dimnames(data) <- list(1:nt, c("X", "Y1", "Y2"))
  if (n[2] > 0) data[(cs[1]+1):cs[2],"Y1"] <- NA
  if (n[3] > 0) data[(cs[2]+1):cs[3],"Y2"] <- NA
  if (n[4] > 0) data[(cs[3]+1):cs[4],c("Y1","Y2")] <- NA
  return(data)
}
# Impute data (from Buuren2018 code in chapter 4)
impute <- function(data, m = 5, method = "norm",
                   print = FALSE, maxit = 10, ...) {
  statistic <- matrix(NA, nrow = maxit, ncol = m)
  for (iter in 1:maxit) {
    if (iter==1) imp <- mice(data, m = m, method = method,
                             print = print, maxit = 1, ...)
    else imp <- mice.mids(imp, maxit = 1, print = print, ...)
    statistic[iter, ] <- unlist(with(imp, cor(Y1, Y2))$analyses)
  }
  return(list(imp = imp, statistic = statistic))
}

# Impute data and save the statistic of interest (estimated correlation between Y1 and Y2)
save_statistics <- function(dt, m = 5, maxit = 10){
  r <- impute(dt, m = m, maxit = maxit)
  rY1Y2 <- c(r$statistic)
  return(rY1Y2)
}

# Simulation --------------------------------------------------------------

# Conditions and specs
  corr_coef <- seq(.94, .99, by = .01)
  m = 5; maxit = 150
  A <- create_cor(corr_coef) # covariance matrices
# Create data
  datasets <- lapply(A, generate, n = c(1000, 4500, 4500, 0)) #90% missings (1e4 cases, 1e3 missings)
    lapply(datasets, md.pattern)
# Perform simulation
  rY1Y2 <- lapply(datasets, save_statistics, m=m, maxit=maxit)

# Plot results ------------------------------------------------------------
# Create result object
  results <- cbind(expand.grid(rep(seq(1, length(A)), each = m*maxit )), #kth conditions
                   rep(rep(seq(1, m), each=maxit), length(A)), # m indicator
                   rep(seq(1, maxit), length(A)*m), #iteration indicator
                   c(do.call(cbind, rY1Y2)) ) # (the statistic of interest)
  colnames(results) <- c("k", "m", "iter", "rY1Y2")

# Create plot
  labels <- paste(as.character(corr_coef), "true Y1-Y2 corr")
  xyplot(rY1Y2 ~ iter | as.factor(k), group = m,
         data = results, layout = c(3,2),
         type="l", as.table = TRUE,
         ylab = "Correlation between Y1 and Y2",
         xlab = "Iteration", col = mdc(3),
         scales = list(y = list(alternating = 1, tck = c(1, 0))),
         strip = strip.custom(bg = "grey95", style = 1,
                              factor.levels = labels))
