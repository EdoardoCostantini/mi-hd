### Title:    Replication Deng Et Al 2016 - Compare methods
### Author:   Edoardo Costantini
### Created:  2020-02-10
### Modified: 2020-02-14

rm(list = ls())

source("./functions.R")
source("./init.R")

out_old <- readRDS("../output/pooled_ZL2016-mc-20200512_1139.rds")
out <- readRDS("../output/pooled_ZL2016-mc-20200512_1145.rds")  # different seed
out <- readRDS("../output/pooled_DEA2016-mc-20200512_2020.rds") # Including Bayesian lasso

# Simulation description
str(out[[1]]$cond_200_4$parms)

# Results from simulation -------------------------------------------------

# Time
res_time <- NULL
for (i in 1:out[[1]]$cond_200_4$parms$dt_rep) {
  res_time <- rbind(res_time, out[[i]]$cond_200_4$run_time_min)
}
round(colMeans(res_time), 3)

# Rsults on a different seed
res_old <- lapply(1:4, function(x){ 
  extract_results(cond_name = names(out_old[[1]])[x], 
                  output = out_old, 
                  dt_rep = 500)
})
names(res_old) <- names(out[[1]]) 
res_old

# New
res <- lapply(1:4, function(x){
  extract_results(cond_name = names(out[[1]])[x], 
                  output = out, 
                  dt_rep = out[[1]]$cond_200_4$parms$dt_rep)
})
names(res) <- names(out[[1]]) 
res

cbind(res_old[[1]], res[[1]])


# Convergence -------------------------------------------------------------
length(out)

# Make a function in the future
dt_rep <- 10
plot(out[[dt_rep]]$cond_200_4$imp_values$MICE_TR)
imp_meth <- c("DURR", "IURR", "blasso")
chain1 <- 1
par(mfrow = c(2, 3))
for (m in 1:length(imp_meth)) {
  for (v in 1:length(parms$z_m_id)) {
    imps_4plot <- out[[dt_rep]]$cond_200_4$imp_values[[imp_meth[m]]][[chain1]][[v]]
    mean_imp <- colMeans(imps_4plot)
    plot(1:parms$iters, mean_imp, type = "l",
         main = paste0(imp_meth[m], " Mean Imputations"),
         ylim = c(mean(mean_imp)-4*sd(mean_imp), mean(mean_imp)+4*sd(mean_imp)),
         ylab = paste0("z", v), xlab = "Iteration")
    for (i in 1:(parms$chains-1)) {
      imps_4plot <- out[[dt_rep]]$cond_200_4$imp_values[[imp_meth[m]]][[i]][[v]]
      mean_imp <- colMeans(imps_4plot)
      lines(1:parms$iters, mean_imp)
    }
  }
}

# Results from paper ------------------------------------------------------
results_paper <- list(cond_200_4 = data.frame(bias_p = c(.074, .017, -.023, 
                                                       -.008, -.356, -.206),
                                              ci_p = c(.894, .930, .946,
                                                     .938, .378, .630)),
                      cond_1000_4 = data.frame(bias_p = c(.066, .030, -.030, 
                                                        -.008, -.356, -.206),
                                               ci_p = c(.904, .918, .934,
                                                      .938, .378, .630)),
                      cond_200_20 = data.frame(bias_p = c(.053, -.012, -.037, 
                                                        -.046, -.264, -.224),
                                               ci_p = c(.886, .918, .942,
                                                      .940, .530, .412)),
                      cond_1000_20 = data.frame(bias_p = c(.032, -.003, -.056, 
                                                         -.046, -.264, -.224),
                                                ci_p = c(.914, .934, .956,
                                                       .940, .530, .412))
)

# Compare paper and simulation study --------------------------------------

compare_b1 <- vector("list", length(out[[1]]))
  names(compare_b1) <- names(out[[1]]) 
for (i in 1:length(out[[1]])) {
  compare_b1[[i]] <- cbind(results_paper[[i]], res_old[[i]], res[[i]])
}
compare_b1
