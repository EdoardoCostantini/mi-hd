### Title:    Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-05-19

rm(list = ls())

source("./functions.R")
source("./init.R")

# Current simulation results ----------------------------------------------
# If you run the simulation script you can use this directly to quick check
# Time
res_time <- NULL
for (i in 1:out[[1]]$cond_200_4$parms$dt_rep) {
  res_time <- rbind(res_time, out[[i]]$cond_200_4$run_time_min)
}
round(colMeans(res_time), 3)

# bias coverage
res <- lapply(1:4, function(x){
  extract_results(cond_name = names(out[[1]])[x], 
                  output = out, 
                  dt_rep = out[[1]]$cond_200_4$parms$dt_rep)
})
names(res) <- names(out[[1]]) 
res

# Convegence
# Make a function in the future
dt_rep <- 1
plot(out[[dt_rep]]$cond_200_4$imp_values$MICE_TR)
imp_meth <- c("DURR", "IURR", "blasso")
chain1 <- 1
par(mfrow = c(3, 3))
for (m in 1:length(imp_meth)) {
  for (v in 1:length(parms$z_m_id)) {
    imps_4plot <- out[[dt_rep]]$cond_200_4$imp_values[[imp_meth[m]]][[chain1]][[v]]
    mean_imp <- rowMeans(imps_4plot)
    plot(1:parms$iters, mean_imp, type = "l",
         main = paste0(imp_meth[m], " Mean Imputations"),
         ylim = c(mean(mean_imp)-4*sd(mean_imp), mean(mean_imp)+4*sd(mean_imp)),
         ylab = paste0("z", v), xlab = "Iteration")
    for (i in 1:(parms$chains-1)) {
      imps_4plot <- out[[dt_rep]]$cond_200_4$imp_values[[imp_meth[m]]][[i]][[v]]
      mean_imp <- rowMeans(imps_4plot)
      lines(1:parms$iters, mean_imp)
    }
  }
}
