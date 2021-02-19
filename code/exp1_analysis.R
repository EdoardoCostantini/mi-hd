### Title:    Analysis of results
### Project:  Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-07-09
### Notes:    reads output form results.R script and shows the numbers that
###           are used to draw the conclusions.

  rm(list = ls())
  source("./init_general.R")

# Use output not saved
  # res <- output
  
# Read results from a run of simulation study
  res <- exp1_res <- readRDS("../output/exp1_simOut_20200731_1735_res.rds")
  res <- exp1_res <- readRDS("../output/exp1_simOut_20200801_1620_res.rds")
  # 1 and 2 are equivalent, but second file has more repetitions (500 vs 750)
  res <- exp1_res <- readRDS("../output/exp1_simOut_20201130_1006_res.rds") # draft 1
  res <- exp1_res <- readRDS("../output/exp1_simOut_2020121516_res.rds") # run with cv bridge (draft 2)
  
# Plot Sizes Decisions

  gp_width <- 15
  gp_height <- 20
  sp_width <- 3.75
  sp_height <- 6

# Conditions
  
  res$conds
  cond_names <- paste0(letters[1:nrow(res$conds)], ") ",
                       "Condition ", 1:nrow(res$conds), ": ",
                       "p = ",  res$conds$p, "; ",
                       "pm = ",  res$conds$pm)
  data.frame(res$conds, cond_names)


# Bias (Facet grid) -------------------------------------------------------
  
  meth_compare = rev(c("DURR_la", "IURR_la", 
                       "blasso", "bridge",
                       "MI_PCA",
                       "MI_CART", "MI_RF", 
                       "missFor", 
                       "CC",
                       "MI_OP"))
  
# > Summary version ####
  pf <- plot_fg(dt = lapply(1:length(res$sem),
                            function(x) data.frame( res$sem[[x]]$bias_per)),
                parPlot = list(means = 1:6,
                               variances = 7:12,
                               covariances = 13:27),
                dt_reps = 1e3, 
                ci_lvl = .95,
                type = "bias",
                summy = TRUE,
                meth_compare = meth_compare)
  pf
  ggsave(file = "../output/graphs/exp1_bias_summy.pdf",
         width = sp_width*4, height = sp_height*3,
         units = "cm",
         pf)
  
# > Supplementary Material ####
  pf <- plot_fg(dt = lapply(1:length(res$sem),
                            function(x) data.frame( res$sem[[x]]$bias_per)),
                parPlot = list(means = 1:6,
                               variances = 7:12,
                               covariances = 13:27),
                dt_reps = 1e3, 
                ci_lvl = .95,
                type = "bias",
                meth_compare = meth_compare)
  pf
  ggsave(file = "../output/graphs/exp1_bias.pdf",
         width = gp_width, height = gp_height,
         units = "cm",
         pf)

# CI (Facet grid) ---------------------------------------------------------
  
  meth_compare = rev(c("DURR_la", "IURR_la", 
                       "blasso", "bridge",
                       "MI_PCA",
                       "MI_CART", "MI_RF", 
                       "missFor", 
                       "CC",
                       "MI_OP"))
# > Summary version ####
  pf <- plot_fg(dt = lapply(1:length(res$sem),
                            function(x) data.frame( res$sem[[x]]$ci_cov)),
                type = "ci",
                parPlot = list(means = 1:6,
                               variances = 7:12,
                               covariances = 13:27),
                dt_reps = 1e3,
                ci_lvl = .95,
                summy = TRUE,
                meth_compare = meth_compare)
  pf
  ggsave(file = "../output/graphs/exp1_CI_summy.pdf",
         width = sp_width*4, height = sp_height*3,
         units = "cm",
         pf)
  
# > Supplementary Material ####
  pf <- plot_fg(dt = lapply(1:length(res$sem),
                            function(x) data.frame( res$sem[[x]]$ci_cov)),
                type = "ci",
                parPlot = list(means = 1:6,
                               variances = 7:12,
                               covariances = 13:27),
                dt_reps = 1e3,
                ci_lvl = .95,
                meth_compare = meth_compare)
  pf
  ggsave(file = "../output/graphs/exp1_CI.pdf",
         width = gp_width, height = gp_height,
         units = "cm",
         pf)
  