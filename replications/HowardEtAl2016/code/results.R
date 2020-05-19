### Title:    Replication Howard et al 2015
### Author:   Edoardo Costantini
### Created:  2020-05-18

rm(list = ls())

source("./functions.R")
source("./init.R")

# Current simulation results ----------------------------------------------

extract_results(cond_name = names(out[[1]])[2], 
                output = out, 
                dt_rep = out[[1]]$cond_0.3_75$parms$dt_rep)

res <- lapply(1:nrow(conds), function(x){
  extract_results(cond_name = names(out[[1]])[x], 
                  output = out, 
                  dt_rep = out[[1]]$cond_0.3_75$parms$dt_rep)
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

# Results from simulation -------------------------------------------------
# Compare results across different simulations (different seeds, set ups etc)

# DURR IURR BLASSO RF TRUE GS CC from Blade
out_blade <- readRDS("../output/pooled_DEA2016-mc-20200517_1023.rds")

# DURR IURR BLASSO TRUE GS CC from Mac
out_bl <- readRDS("../output/pooled_DEA2016-mc-20200513_1532.rds") # old gibbs shape
# out_bl <- readRDS("../output/pooled_DEA2016-mc-20200516_1712.rds") # mac run after reshape

# ZahoLong replication compariosn
out_ZL2016 <- readRDS("../../ZhaoLong2016/output/pooled_ZL2016-mc-20200507_1521.rds") # run with ZL2016 code

  
## --------- ##
## Blade Run ##
## --------- ## 
  
# Time
  out_blade_time <- NULL
  for (i in 1:out_blade[[1]]$cond_200_4$parms$dt_rep) {
    out_blade_time <- rbind(out_blade_time, out_blade[[i]]$cond_200_4$run_time_min)
  }
  
# Bias / coverage
  out_blade_res <- extract_results(cond_name = names(out_blade[[1]]), 
                                   output = out_blade, 
                                   dt_rep = 500)
  names(out_blade[[1]]) # single condition
  
  out_blade_res
  data.frame(minutes = round(colMeans(out_blade_time), 3))
  
# Imputation Convergence
  dt_rep <- 5 # for which dataset
  plot_iters <- out_blade[[dt_rep]]$cond_200_4$parms$iters
  plot_chain <- out_blade[[dt_rep]]$cond_200_4$parms$chains
  plot_misvar <- out_blade[[dt_rep]]$cond_200_4$parms$z_m_id
  imp_meth <- c("DURR", "IURR", "blasso")
  chain1 <- 1
  par(mfcol = c(3, 3))
  for (m in 1:length(imp_meth)) {
    for (v in 1:length(plot_misvar)) {
      imps_4plot <- out_blade[[dt_rep]]$cond_200_4$imp_values[[imp_meth[m]]][[chain1]][[v]]
      mean_imp <- rowMeans(imps_4plot)
      plot(1:plot_iters, mean_imp, type = "l",
           main = paste0(imp_meth[m], " Mean Imputations"),
           ylim = c(mean(mean_imp)-4*sd(mean_imp), mean(mean_imp)+4*sd(mean_imp)),
           ylab = paste0("z", v), xlab = "Iteration")
      for (i in 1:(plot_chain-1)) {
        imps_4plot <- out_blade[[dt_rep]]$cond_200_4$imp_values[[imp_meth[m]]][[i]][[v]]
        mean_imp <- rowMeans(imps_4plot)
        lines(1:plot_iters, mean_imp)
      }
    }
  }
  plot(out_blade[[dt_rep]]$cond_200_4$imp_values$MICE_TR)

## --------- ##
##  Mac Run  ##
## --------- ## 
## Note: Pre reshaping of samplers
  
# Time
  out_bl_time <- NULL
  for (i in 1:out_bl[[1]]$cond_200_4$parms$dt_rep) {
    out_bl_time <- rbind(out_bl_time, out_bl[[i]]$cond_200_4$run_time_min)
  }
  
# Bias / coverage
  out_bl_res <- extract_results(cond_name = names(out_bl[[1]]), 
                                output = out_bl, 
                                dt_rep = 200)
  names(out_bl[[1]]) # single condition
  
  out_bl_res
  round(colMeans(out_bl_time), 3)
  
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

## --------- ##
## Univa Run ##
## --------- ## 
## Note: For *all* analysis model parameters
  out_ZL2016_res <- lapply(1:4, function(x){
    extract_results(cond_name = names(out_ZL2016[[1]])[x], 
                    output = out_ZL2016, 
                    dt_rep = out_ZL2016[[1]]$cond_200_4$parms$dt_rep)
  })
  names(out_ZL2016_res) <- names(out_ZL2016[[1]]) 
  out_ZL2016_res
