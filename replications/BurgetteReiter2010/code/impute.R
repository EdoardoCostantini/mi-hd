### Title:    Replication Burgette Reiter 2010 - Perform imputation loops
### Author:   Edoardo Costantini
### Created:  2020-JAN-10
### Modified: 2020-JAN-10
### Descript: given 1e3 datasets stored in one large matrix, this script
###           imputes the values according to seqeuntial CART MICE and
###           bayesian linear regresion MICE.
### input: many datasets w/ missing values stored in a large matrix
###        of dimensions (n * rep_datasets) x (p * 2)
### output: a list containing (1) estimates of regression coefficeints pooled 
###         from m imputed datasest, (2) pooled confidence intervals, (3) out
###         of sample rmse

    library(mvtnorm) # mvnorm
    library(mice)    # ampute
    library(mi)
    
    source("./functions.R")                  # from source file location
    source("../../../R/functions_allpurp.R") # for bbsample function
    
# Set up ------------------------------------------------------------------

  # Get data
  dt_matrix <- readRDS("../output/data_2020-02-13_rep_50.rds")
  p <- ncol(dt_matrix)/2
  
  # Iterations  
  iter <- 1e3
  output_List <- vector("list", 3)
  
for (i in 1:iter) {
  
  time.start <- Sys.time()

# Select data -------------------------------------------------------------
  
  row_indx <- ((i-1)*1e3 + 1):( (i)*1e3 ) 
    # where in the matrix is the i-th dataset
  
  Xy_mis <- dt_matrix[row_indx, 1:p]
  Xy_test <- dt_matrix[row_indx, (p+1):(2*p)]
  X_test <- Xy_test[,-p]
  y_test <- Xy_test[,p]
  
# Imputation --------------------------------------------------------------

  # MICE cart Burgette Reiter
  imp_CART_bb <- miceImpHDv::parlmice(Xy_mis, 
                                      # Parallel details
                                      n.core = 10,
                                      n.imp.core = 1, # in each core, perform 1 
                                                      # imputation of m desired 
                                                      # ones
                                      cl.type = "FORK",
                                      cluster.seed = 20200211,
                                      
                                      # Iteration details
                                      #m = n_cores * n.imp.core
                                      maxit = 10, # l in the paper
                                      
                                      # Method details
                                      meth = 'cart.bb', 
                                      minbucket = 5) # minimum leaf size

  # MICE default
  imp_norm <- mice::parlmice(Xy_mis,
                                n.core = 10,
                                n.imp.core = 1, 
                                cl.type = "FORK",
                                cluster.seed = 20200211,
                                maxit = 10,
                                meth = 'norm') # norm is the deafult method
                                               # B&R 2010 are referring to
  
# Pooling -----------------------------------------------------------------
  
  # MICE cart Burgette Reiter
  fit <- with(imp_CART_bb, 
              expr = lm(V11 ~ V1 + V2 + V3 + V8 + V9 + I(V3^2) + V1:V2 + V8:V9))
                      # V11 = y
  pool_CART_est <- mice::pool(fit)$pooled[,1]
  
  pool_CART_conf <- summary(mice::pool(fit), 
                            conf.int = TRUE)[,c("2.5 %", "97.5 %")]
  
  # MICE default
  fit_norm <- with(imp_norm, 
              expr = lm(V11 ~ V1 + V2 + V3 + V8 + V9 + I(V3^2) + V1:V2 + V8:V9))
  pool_norm_est <- mice::pool(fit_norm)$pooled[,1]
  
  pool_norm_conf <- summary(mice::pool(fit_norm), 
                            conf.int = TRUE)[,c("2.5 %", "97.5 %")]
  
  # group estiamtes
  pool_est <- cbind(pool_CART_est, pool_norm_est)
  pool_conf <- list(pool_CART_conf = pool_CART_conf, 
                    pool_mi_conf = pool_norm_conf)

# Prediction --------------------------------------------------------------
  
  y_hat <- apply(pool_est, 2, gen_yhat_BR1073, x = X_test)
    colnames(y_hat) <- c("CART", "mi")

# Save output -------------------------------------------------------------

  output <- list(pool_est = pool_est,
                 pool_conf = pool_conf,
                 y_hat = y_hat,
                 y_test = y_test)
  output_List[[1]][[i]] <- output$pool_est
  output_List[[2]][[i]] <- data.frame(output$y_hat,
                                      output$y_test)
  output_List[[3]][[i]] <- output$pool_conf
  
  # Monitor progress
  print(paste0("Repetition ", i, " took: ", round(Sys.time() - time.start, 1)))
  
}
  
dirOut <- "../output/"
fileName <- paste0("imp_out_", Sys.time(), "_rep_", iter, ".rds")

saveRDS(output_List, paste0(dirOut, fileName))
  
  