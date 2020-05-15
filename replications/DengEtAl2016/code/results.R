### Title:    Replication Deng Et Al 2016 - Compare methods
### Author:   Edoardo Costantini
### Created:  2020-02-10
### Modified: 2020-02-14

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

# Comparison of different simulaations ------------------------------------
# Compare results across different simulations (different seeds, set ups etc)

# DURR IURR exclusively (replication ZhaoLong2016)
out_s1 <- readRDS("../output/pooled_ZL2016-mc-20200512_1139.rds") # seed 1
out_s2 <- readRDS("../output/pooled_ZL2016-mc-20200512_1145.rds") # seed 2
  # Info
  out_s1[[1]]$cond_200_4$parms$detlamod
  out_s2[[1]]$cond_200_4$parms$detlamod
  out_s2[[1]]$cond_200_4$cond_bias

# DURR IURR BLASSO
# out <- readRDS("../output/pooled_DEA2016-mc-20200512_2020.rds") # Including Bayesian lasso
out_bl <- readRDS("../output/pooled_DEA2016-mc-20200513_1532.rds") # run again
  # Info
  out_bl[[1]]$cond_200_4$parms$detlamod
  out_bl[[1]]$cond_200_4$cond_bias

# Random Forest
out_rf <- readRDS("../output/pooled_DEA2016-mc-20200513_1211.rds") # Just Rf 50 reps
  # Info
  out_rf[[1]]$cond_200_4$parms$dt_rep
  out_rf[[1]]$cond_200_4$parms$detlamod
  out_rf[[1]]$cond_200_4$cond_bias

# ZahoLong replication compariosn
out_ZL2016 <- readRDS("../../ZhaoLong2016/output/pooled_ZL2016-mc-20200507_1521.rds") # run with ZL2016 code
out_ZL2016_update <- readRDS("../output/pooled_ZL2016-mc-20200513_1833.rds") # run with DEAT2016 code
  # Info
  out_ZL2016[[1]]$cond_200_4$parms$dt_rep
  out_ZL2016_update[[1]]$cond_200_4$parms$detlamod
  out_ZL2016_update[[1]]$cond_200_4$cond_bias 
    # naming issue: MI_T is MICE_RF and move everything else

# Results from simulation -------------------------------------------------

# 1. Multivairate check DURR IURR blasso ####
  
# Time
out_bl_time <- NULL
for (i in 1:out_bl[[1]]$cond_200_4$parms$dt_rep) {
  out_bl_time <- rbind(out_bl_time, out_bl[[i]]$cond_200_4$run_time_min)
}
round(colMeans(out_bl_time), 3)

# Bias / coverage
out_bl_res <- extract_results(cond_name = names(out_bl[[1]]), 
                              output = out_bl, 
                              dt_rep = 200)
names(out_bl[[1]]) # single condition
out_bl_res

# Imputation Convergence
dt_rep <- 3 # for which dataset

plot_iters <- out_bl[[dt_rep]]$cond_200_4$parms$iters
plot_chain <- out_bl[[dt_rep]]$cond_200_4$parms$chains
plot_misvar <- out_bl[[dt_rep]]$cond_200_4$parms$z_m_id
imp_meth <- c("DURR", "IURR", "blasso")
chain1 <- 1
par(mfcol = c(3, 3))
for (m in 1:length(imp_meth)) {
  for (v in 1:length(plot_misvar)) {
    imps_4plot <- out_bl[[dt_rep]]$cond_200_4$imp_values[[imp_meth[m]]][[chain1]][[v]]
    mean_imp <- colMeans(imps_4plot)
    plot(1:plot_iters, mean_imp, type = "l",
         main = paste0(imp_meth[m], " Mean Imputations"),
         ylim = c(mean(mean_imp)-4*sd(mean_imp), mean(mean_imp)+4*sd(mean_imp)),
         ylab = paste0("z", v), xlab = "Iteration")
    for (i in 1:(plot_chain-1)) {
      imps_4plot <- out_bl[[dt_rep]]$cond_200_4$imp_values[[imp_meth[m]]][[i]][[v]]
      mean_imp <- colMeans(imps_4plot)
      lines(1:plot_iters, mean_imp)
    }
  }
}
plot(out_bl[[dt_rep]]$cond_200_4$imp_values$MICE_TR)

# 2. Multivairate check Random Forest ####

# Time
out_rf_time <- NULL
for (i in 1:out_rf[[1]]$cond_200_4$parms$dt_rep) {
  out_rf_time <- rbind(out_rf_time, out_rf[[i]]$cond_200_4$run_time_min)
}
round(colMeans(out_rf_time), 3)

# Bias / coverage
out_rf_res <- extract_results(cond_name = names(out_rf[[1]]), 
                              output = out_rf, 
                              dt_rep = 50)
names(out_rf[[1]]) # single condition
out_rf_res
out_bl_res$bias
out_rf_res$bias

# Imputation Convergence
dt_rep <- 3 # for which dataset
plot(out_rf[[dt_rep]]$cond_200_4$imp_values$MICE_TR)
plot(out_rf[[dt_rep]]$cond_200_4$imp_values$MICE_RF) #  does not work: change object saved

# 3. Univariate comparison ####
# Equivalence of ZhaoLong2016 replication with DEAT2016 code

# bias coverage
out_ZL2016_res <- lapply(1:2, function(x){
  extract_results(cond_name = names(out_ZL2016[[1]])[x], 
                  output = out_ZL2016, 
                  dt_rep = out_ZL2016[[1]]$cond_200_4$parms$dt_rep)
})
names(out_ZL2016_res) <- names(out_ZL2016_update[[1]]) 
out_ZL2016_res

out_ZL2016_up_res <- lapply(1:2, function(x){
  extract_results(cond_name = names(out_ZL2016_update[[1]])[x], 
                  output = out_ZL2016_update, 
                  dt_rep = out_ZL2016_update[[1]]$cond_200_4$parms$dt_rep)
})
names(out_ZL2016_up_res) <- names(out_ZL2016_update[[1]]) 
out_ZL2016_up_res

comp_oldnew <- data.frame(
  old = out_ZL2016_res$cond_200_4$bias[c("DURR", "IURR", "blasso", "MI_T", "CC"), "z1"],
  new = out_ZL2016_up_res$cond_200_4$bias[c(1, 2, 3, 5, 7),
                                    "z1"]
)

row.names(comp_oldnew) <- c("DURR", "IURR", "blasso", "MI_T", "CC")
comp_oldnew
