### Title:    Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-05-19

source("./init.R")
source("./functions.R")

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
