### Title:    Analysis of results from experiment 2
### Project:  Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-07-09
### Notes:    reads output form results.R script and shows the numbers that
###           are used to draw the conclusions.

  rm(list = ls())
  source("./init_general.R")
  
  # Read results from a run of simulation study
  # filename <- "exp2_simOut_20201116_1743" # short run with scaslde items
  # filename <- "exp2_simOut_20201117_1125" # short run with scaslde items
  # filename <- "exp2_simOut_20200819_1743" # last good one

  # LVs as MAR (Joined file form two sources) 1e3 reps
  filename <- "exp2_simOut_2020122417"
  res <- readRDS(paste0("../output/", filename, "_res.rds"))

  # Items as MAR 1e3 reps
  filename <- "exp2_simOut_20210728_1351_res"
  res <- readRDS(paste0("../output/", filename, ".rds"))

  # Plot Sizes Decisions
  gp_width <- 15
  gp_height <- 20
  sp_width <- 3.75
  sp_height <- 6
  
# Conds
  res$conds

  # Full information names
  cond_parms <- paste0("l = ",  res$conds$lv, ", ",
                       "pm = ",  res$conds$pm, "\n",
                       "\u03bb = ", c(rep("(0.9, 0.97)", 4),
                                      rep("(0.5, 0.6)", 4)))
  
  cond_labels <- c("low-dim-low-pm-high-\u03bb",
                  "high-dim-low-pm-high-\u03bb",
                  "low-dim-high-pm-high-\u03bb",
                  "high-dim-high-pm-high-\u03bb",
                  "low-dim-low-pm-low-\u03bb",
                  "high-dim-low-pm-low-\u03bb",
                  "low-dim-high-pm-low-\u03bb",
                  "high-dim-high-pm-low-\u03bb")
  cond_names <- paste(cond_labels, cond_parms, sep = " \n ")
  # Numbered Names
  # cond_names <- paste0("Condition ", 1:8)
  data.frame(res$conds, cond_names)
  
# Select conditions to print
  conds_select <- 1:4
  conds_select <- 5:8

# Bias (Facet grid) -------------------------------------------------------

  meth_compare = rev(c("DURR_la", "IURR_la", 
                       "blasso", "bridge",
                       "MI_PCA",
                       "MI_CART", "MI_RF", 
                       "missFor", 
                       "CC",
                       "MI_OP"))
  
# > Summary version ####
  pf <- plot_fg(dt = lapply(1:length(res$semR),
                            function(x) data.frame( res$semR[[x]]$bias_per))[conds_select],
                type = "bias",
                parPlot = list(means = 1:10,
                               variances = 11:20,
                               covariances = 21:65),
                dt_reps = 1e3,
                ci_lvl = .95,
                cond_labels = cond_names[conds_select],
                summy = TRUE,
                meth_compare = meth_compare)
  pf
  # Save Plot
  ggsave(file = paste0("../output/graphs/",
                       "exp2_semR_bias_",
                       paste0(range(conds_select), collapse = ""),
                       "_summy",
                       ".pdf"),
         width = sp_width*4, height = sp_height*3,
         units = "cm",
         device = cairo_pdf, # allows for greek letters in text
         pf)
  
# > Supplementary Material ####
  pf <- plot_fg(dt = lapply(1:length(res$semR),
                            function(x) data.frame( res$semR[[x]]$bias_per))[conds_select],
                type = "bias",
                parPlot = list(means = 1:10,
                               variances = 11:20,
                               covariances = 21:65),
                dt_reps = 1e3,
                ci_lvl = .95,
                cond_labels = cond_names[conds_select],
                summy = FALSE,
                meth_compare = meth_compare)
  pf
  
  # Save Plot
  ggsave(file = paste0("../output/graphs/",
                       "exp2_semR_bias_",
                       paste0(range(conds_select), collapse = ""),
                       ".pdf"),
         width = gp_width, height = gp_height,
         units = "cm",
         device = cairo_pdf, # allows for greek letters in text
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
  pf <- plot_fg(dt = lapply(1:length(res$semR),
                            function(x) data.frame( res$semR[[x]]$ci_cov))[conds_select],
                type = "ci",
                parPlot = list(means = 1:10,
                               variances = 11:20,
                               covariances = 21:65),
                dt_reps = 1e3,
                ci_lvl = .95,
                cond_labels = cond_names[conds_select],
                summy = TRUE,
                meth_compare = meth_compare)
  pf
  ggsave(file = paste0("../output/graphs/",
                       "exp2_semR_ci_",
                       paste0(range(conds_select), collapse = ""),
                       "_summy",
                       ".pdf"),
         width = sp_width*4, height = sp_height*3,
         units = "cm",
         device = cairo_pdf, # allows for greek letters in text
         pf)
# > Supplementary Material ####
  pf <- plot_fg(dt = lapply(1:length(res$semR),
                            function(x) data.frame( res$semR[[x]]$ci_cov))[conds_select],
                type = "ci",
                parPlot = list(means = 1:10,
                               variances = 11:20,
                               covariances = 21:65),
                dt_reps = 1e3,
                ci_lvl = .95,
                cond_labels = cond_names[conds_select],
                meth_compare = meth_compare)
  pf
  ggsave(file = paste0("../output/graphs/",
                       "exp2_semR_ci_",
                       paste0(range(conds_select), collapse = ""),
                       ".pdf"),
         width = gp_width, height = gp_height,
         units = "cm",
         device = cairo_pdf, # allows for greek letters in text
         pf)

# CFA (Facet grid) --------------------------------------------------------
  
  meth_compare = rev(c("DURR_la", "IURR_la", 
                       "blasso", "bridge",
                       "MI_PCA",
                       "MI_CART", "MI_RF", 
                       "missFor", 
                       "CC",
                       "MI_OP"))
  
  # > Summary version ####
  # Conditona 1-4
  pt.1 <- plot_fg(dt = lapply(1:length(res$CFA),
                              function(x) data.frame( res$CFA[[x]]$bias_per))[1:4],
                  type = "bias",
                  parPlot = list(Loadings = 1:10),
                  dt_reps = 1e3,
                  ci_lvl = .95,
                  cond_labels = cond_names[1:4],
                  summy = TRUE,
                  meth_compare = meth_compare)
  
  # Conditions 5-8
  pt.2 <- plot_fg(dt = lapply(1:length(res$CFA),
                            function(x) data.frame( res$CFA[[x]]$bias_per))[5:8],
                type = "bias",
                parPlot = list(Loadings = 1:10),
                dt_reps = 1e3,
                ci_lvl = .95,
                cond_labels = cond_names[5:8],
                summy = TRUE,
                meth_compare = meth_compare)
  
  # Combine Elements
  pf <- cowplot::ggdraw() +
    cowplot::draw_plot(pt.1, 
                       x = 0, y = .5, 
                       width = 1, height = .5) +
    cowplot::draw_plot(pt.2, 
                       x = 0, y = 0, 
                       width = 1, height = .5) +
    cowplot::draw_plot_label(label = c("(A)", "(B)"), 
                             x = c(0, 0), 
                             y = c(1, .5), 
                             size = 5,
                             fontface = "bold")
  pf
  
  # Save Plot
  ggsave(file = "../output/graphs/exp2_CFA_lambda_BPR_summy.pdf",
         width = sp_width*4, height = sp_height*2,
         units = "cm",
         device = cairo_pdf, # allows for greek letters in text
         pf)

  # > Supplementary Material: Bias ####
  # Conditona 1-4
  pt.1 <- plot_fg(dt = lapply(1:length(res$CFA),
                              function(x) data.frame( res$CFA[[x]]$bias_per))[1:4],
                  type = "bias",
                  parPlot = list(Loadings = 1:10),
                  dt_reps = 1e3,
                  ci_lvl = .95,
                  cond_labels = cond_names[1:4],
                  meth_compare = meth_compare)
  
  # Conditions 5-8
  pt.2 <- plot_fg(dt = lapply(1:length(res$CFA),
                              function(x) data.frame( res$CFA[[x]]$bias_per))[5:8],
                  type = "bias",
                  parPlot = list(Loadings = 1:10),
                  dt_reps = 1e3,
                  ci_lvl = .95,
                  cond_labels = cond_names[5:8],
                  meth_compare = meth_compare)
  
  # Combine Elements
  pf <- cowplot::ggdraw() +
    cowplot::draw_plot(pt.1, 
                       x = 0, y = .5, 
                       width = 1, height = .5) +
    cowplot::draw_plot(pt.2, 
                       x = 0, y = 0, 
                       width = 1, height = .5) +
    cowplot::draw_plot_label(label = c("(A)", "(B)"), 
                             x = c(0, 0), 
                             y = c(1, .5), 
                             size = 5,
                             fontface = "bold")
  pf
  
  # Save Plot
  ggsave(file = "../output/graphs/exp2_CFA_lambda_BPR.pdf",
         width = gp_width, height = gp_height*2/3,
         units = "cm",
         device = cairo_pdf, # allows for greek letters in text
         pf)
  
  # > Supplementary Material: CIC ####
  
  # Conditona 1-4
  pt.1 <- plot_fg(dt = lapply(1:length(res$CFA),
                              function(x) data.frame( res$CFA[[x]]$ci_cov))[1:4],
                  type = "ci",
                  parPlot = list(Loadings = 1:10),
                  dt_reps = 1e3,
                  ci_lvl = .95,
                  cond_labels = cond_names[1:4],
                  meth_compare = meth_compare)
  
  # Conditions 5-8
  pt.2 <- plot_fg(dt = lapply(1:length(res$CFA),
                              function(x) data.frame( res$CFA[[x]]$ci_cov))[5:8],
                  type = "ci",
                  parPlot = list(Loadings = 1:10),
                  dt_reps = 1e3,
                  ci_lvl = .95,
                  cond_labels = cond_names[5:8],
                  meth_compare = meth_compare)
  
  # Combine Elements
  pf <- cowplot::ggdraw() +
    cowplot::draw_plot(pt.1, 
                       x = 0, y = .5, 
                       width = 1, height = .5) +
    cowplot::draw_plot(pt.2, 
                       x = 0, y = 0, 
                       width = 1, height = .5) +
    cowplot::draw_plot_label(label = c("(A)", "(B)"), 
                             x = c(0, 0), 
                             y = c(1, .5), 
                             size = 5,
                             fontface = "bold")
  pf
  
  # Save Plot
  ggsave(file = "../output/graphs/exp2_CFA_lambda_CIC.pdf",
         width = gp_width, height = gp_height*2/3,
         units = "cm",
         device = cairo_pdf, # allows for greek letters in text
         pf)