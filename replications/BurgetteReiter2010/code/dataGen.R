### Title:    Replication Burgette Reiter 2010 - Create data
### Author:   Edoardo Costantini
### Created:  2020-JAN-10
### Modified: 2020-JAN-10
### Notes:    Create 1e3 simulated datasets for sim study p.1072

library(mvtnorm) # mvnorm
library(mice)    # ampute
library(mi)

source("./functions.R") # from source file location

# Set up ------------------------------------------------------------------

dirOut <- "../output/"
fileName <- paste0("data_", Sys.Date(), "_rep_", iter, ".rds")

iter <- 1e3 # 1e3 goal datasets

mypatterns <- rbind(c(0,0,0,0,0,0,0,0,1,1,0), # 1 = observed variable
                    c(1,0,0,0,0,0,0,0,1,1,0),
                    c(0,1,0,0,0,0,0,0,1,1,0),
                    c(0,0,1,0,0,0,0,0,1,1,0),
                    c(0,0,0,1,0,0,0,0,1,1,0),
                    c(0,0,0,0,1,0,0,0,1,1,0),
                    c(0,0,0,0,0,1,0,0,1,1,0),
                    c(0,0,0,0,0,0,1,0,1,1,0),
                    c(0,0,0,0,0,0,0,1,1,1,0),
                    c(0,0,0,0,0,0,0,0,1,1,1),
                    c(0,0,1,1,0,0,1,1,1,1,0),
                    c(0,0,0,0,1,1,0,0,1,1,1),
                    c(0,1,0,1,0,1,0,1,1,1,0),
                    c(1,0,1,0,1,0,1,0,1,1,1)
)

p <- ncol(mypatterns)

# Loop --------------------------------------------------------------------

storeData <- matrix(rep(NA, iter*1e3*2*p), ncol = 2*p)

set.seed(20200213)

for (i in 1:iter) {
  # Create imputation data according to:
  # lm(y ~ V1 + V2 + V3 + V8 + V9 + I(V3^2) + V1:V2 + V8:V9, data = Xy)
  X <- gen_X_BR1073(n = 1e3)
  y <- gen_y_BR1073(X)
  Xy <- as.data.frame(cbind(X, y))
  
  # Impose missingness with patterns defined at beginning
  Xy_amp <- ampute(Xy,
                   prop = .75, # percentage of cases w/ missing values
                   mech = "MAR",
                   patterns = mypatterns) 
  Xy_mis <- Xy_amp$amp

  # Create test dataset (fully observed)
  X_test <- gen_X_BR1073(n = 500)
  y_test <- gen_y_BR1073(X_test)
  Xy_test <- as.data.frame(cbind(X_test, y_test))
  
  # Store datasets in a matrix
  row_indx <- ((i-1)*1e3 + 1):( (i)*1e3 ) # where to store in the matrix
  
  storeData[row_indx, 1:p] <- as.matrix(Xy_mis)
  storeData[row_indx, (p+1):(2*p)] <- as.matrix(Xy_test)

  # Print status
  print(paste0("Dataset ", i, " done, ", Sys.time()))
}

saveRDS(storeData, paste0(dirOut, fileName))

