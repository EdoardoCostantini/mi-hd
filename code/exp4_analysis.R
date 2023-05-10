# Title:    Analysis of results from experiment 4 (EVS)
# Porject:  Imputing High Dimensional Data
# Author:   Edoardo Costantini
# Created:  2020-08-27
# Modified: 2023-03-30
# Notes:    reads output form results.R script and shows the numbers that
#           are used to draw the conclusions.

  rm(list = ls())
  source("./init_general.R")

# Read results from a run of simulation study
  # Result selection
  filename <- "exp4_simOut_20201016_2341" # old
  filename <- "exp4_simOut_20201019_1344" # current one
  filename <- "exp4_simOut_20201027_1610" # combined results from more runs
  
  # Current runs
  filename <- "exp4_simOut_20201204_2121" # updated model 1, 500 data
  filename <- "exp4_simOut_20201207_1134" # same seed as 20201204_2121, but next 500 samples
  filename <- "exp4_simOut_2020120704" # joined 20201204_2121 and 20201207_1134
  filename <- "exp4_simOut_20220226_0950" # see readme for details
  filename <- "exp4_simOut_20220131_1603" # joined 20201204_2121 and 20201207_1134 + MI-qp and MI-am
  filename <- "exp4_simOut_20220226_0950" # see readme for details
  filename <- "exp4_simOut_20230323_1551" # see readme for details
  
  # Read R object
  res <- readRDS(paste0("../output/", filename, "_res.rds"))
  out <- res
  
  # Plot Size
  sp_width <- 5
  sp_height <- 4

# Summary of set up -------------------------------------------------------

  # Bias percent
  # M1
  lapply(1:length(res$m1),
         function(x) res$m1[[x]]$bias_per["rel", ])
  # M2
  lapply(1:length(res$m2),
         function(x) res$m2[[x]]$bias_per["NatAt", ])
  
  # Bias sd
  # M1
  lapply(1:length(res$m1),
         function(x) round(res$m1[[x]]$bias_sd, 2)["rel", ] )
  # M2
  lapply(1:length(res$m2),
         function(x) round(res$m2[[x]]$bias_sd, 2)["NatAt", ] )
  
  # Confidence Intervals
  lapply(1:length(res$m1),
         function(x) res$m1[[x]]$ci_cov["rel", ])[[1]]
  lapply(1:length(res$m2),
         function(x) res$m2[[x]]$ci_cov["NatAt", ])
  lapply(1:length(res$m2),
           function(x) res$m1[[x]]$ci_cov)[[1]]

# Variable of interest ----------------------------------------------------

  # > Bias (Facet) ####
  # Which Methods do you want to plot
  meths = rev(c("DURR_la", "IURR_la", 
                "blasso", "bridge",
                "MI_PCA",
                "MI_CART", 
                "MI_RF", 
                "stepFor", 
                "missFor", 
                "CC",
                "MI_qp",
                "MI_am",
                "MI_OP"))
  pf <- plot_exp4(dt = list(lapply(1:length(res$m1),
                                   function(x) res$m1[[x]]$bias_per["rel", ]),
                            lapply(1:length(res$m2),
                                   function(x) res$m2[[x]]$bias_per["NatAt", ])),
            type = "bias",
            dt_reps = 1e3,
            ci_lvl = .95,
            plot_cond = NULL,
            plot_name = NULL,
            bar_col = "darkgray",
            meth_compare = meths,
            meth_sort = FALSE)
  pf
  
  ggsave(pf,
         file = "../output/graphs/exp4_imp_bias.pdf",
         width = sp_width*2, height = sp_height*2,
         units = "cm",
         device = cairo_pdf)

  # > CI (Facet) ####
  # Which Methods do you want to plot
  meths = rev(c("DURR_la", "IURR_la", 
                "blasso", "bridge",
                "MI_PCA",
                "MI_CART", 
                "MI_RF", 
                "stepFor", 
                "missFor", 
                "CC",
                "MI_qp",
                "MI_am",
                "MI_OP",
                "GS"))
  
  pf <- plot_exp4(dt = list(lapply(1:length(res$m1),
                             function(x) res$m1[[x]]$ci_cov["rel", ]),
                      lapply(1:length(res$m2),
                             function(x) res$m2[[x]]$ci_cov["NatAt", ])),
            type = "ci",
            dt_reps = 1e3,
            ci_lvl = .95,
            plot_cond = NULL,
            plot_name = NULL,
            bar_col = "darkgray",
            meth_compare = meths,
            meth_sort = FALSE)
  pf
  ggsave(pf, 
         file = "../output/graphs/exp4_imp_ci.pdf", 
         width = sp_width*2, height = sp_height*2,
         units = "cm",
         device = cairo_pdf)
  
  # > CIW (Facet) ####
  # Which Methods do you want to plot
  meths = rev(c("DURR_la", "IURR_la", 
                "blasso", "bridge",
                "MI_PCA",
                "MI_CART", 
                "MI_RF", 
                "stepFor", 
                "missFor", 
                "CC",
                "MI_qp",
                "MI_am",
                "MI_OP",
                "GS"))
  
  pf <- plot_exp4(dt = list(lapply(1:length(res$m1),
                             function(x) data.frame(t(colMeans(res$m1[[x]]$CIW)))),
                      lapply(1:length(res$m2),
                             function(x) data.frame(t(colMeans(res$m2[[x]]$CIW))))),
            type = "ciw",
            dt_reps = 500,
            ci_lvl = .95,
            plot_cond = NULL,
            plot_name = NULL,
            bar_col = "#595959",
            meth_compare = meths,
            meth_sort = FALSE)
  ggsave(pf, 
         file = "../output/graphs/exp4_imp_ciw.pdf", 
         width = sp_width*2, height = sp_height*2,
         units = "cm",
         device = cairo_pdf)
  
# All parameters ----------------------------------------------------------

  # Which Methods do you want to plot
  meths = rev(c("DURR_la", "IURR_la", 
                "blasso", "bridge",
                "MI_PCA",
                "MI_CART", 
                "MI_RF",
                "stepFor",
                "MI_qp",
                "MI_am",
                "MI_OP",
                "CC"))

  # > Bias (Facet) ####
  pt1 <- plot_exp4_meth(dt = lapply(1:length(res$m1),
                                    function(x) res$m1[[x]]$bias_per),
                        type = "bias",
                        dt_reps = 1e3,
                        ci_lvl = .95,
                        meth_compare = meths)

  pt2 <- plot_exp4_meth(dt = lapply(1:length(res$m2),
                                    function(x) res$m2[[x]]$bias_per),
                        type = "bias",
                        dt_reps = 1e3,
                        ci_lvl = .95,
                        meth_compare = meths)
  
  # Save
  ggsave(pt1,
         file = "../output/graphs/exp4_imp_bias_allParms_m1.pdf",
         width = 15, height = 19,
         units = "cm",
         device = cairo_pdf)
  ggsave(pt2,
         file = "../output/graphs/exp4_imp_bias_allParms_m2.pdf",
         width = 15, height = 19,
         units = "cm",
         device = cairo_pdf)

  # > CI (Facet) ####
  # Which Methods do you want to plot
  meths = rev(c("DURR_la", "IURR_la",
                "blasso", "bridge",
                "MI_PCA",
                "MI_CART",
                "MI_RF",
                "stepFor",
                "MI_qp",
                "MI_am",
                "MI_OP",
                "CC",
                "GS"))

  pt1 <- plot_exp4_meth(dt = lapply(1:length(res$m1),
                                    function(x) res$m1[[x]]$ci_cov),
                        type = "ci",
                        dt_reps = 1e3,
                        ci_lvl = .95,
                        meth_compare = meths)

  pt2 <- plot_exp4_meth(dt = lapply(1:length(res$m2),
                                    function(x) res$m2[[x]]$ci_cov),
                        type = "ci",
                        dt_reps = 1e3,
                        ci_lvl = .95,
                        meth_compare = meths)

  ggsave(pt1,
         file = "../output/graphs/exp4_imp_ci_allParms_m1.pdf",
         width = 15, height = 21,
         units = "cm",
         device = cairo_pdf)
  ggsave(pt2,
         file = "../output/graphs/exp4_imp_ci_allParms_m2.pdf",
         width = 15, height = 21,
         units = "cm",
         device = cairo_pdf)

# > CIW (Facet) ####
pt1 <- plot_exp4_meth(
       dt = lapply(
              1:length(res$m1),
              function(x) res$m1[[x]]$CIW
       ),
       type = "CIW",
       dt_reps = 1e3,
       ci_lvl = .95,
       meth_compare = meths
)
pt2 <- plot_exp4_meth(
       dt = lapply(
              1:length(res$m2),
              function(x) res$m2[[x]]$CIW
       ),
       type = "CIW",
       dt_reps = 1e3,
       ci_lvl = .95,
       meth_compare = meths
)

  ggsave(pt1,
         file = "../output/graphs/exp4_imp_ciw_allParms_m1.pdf",
         width = 15, height = 21,
         units = "cm",
         device = cairo_pdf)
  ggsave(pt2,
         file = "../output/graphs/exp4_imp_ciw_allParms_m2.pdf",
         width = 15, height = 21,
         units = "cm",
         device = cairo_pdf)
  
# Multivariate distance ---------------------------------------------------
  
  # > ED BIAS (Facet) ####
  # Which Methods do you want to plot
  meths = rev(c("DURR_la", "IURR_la", 
                "blasso", "bridge",
                "MI_PCA",
                "MI_CART", "MI_RF", 
                "missFor", 
                "CC",
                "MI_OP"))
  
  pf <- plot_exp4(dt = list(lapply(1:length(res$m1),
                             function(x) res$m1[[x]]$ed_est),
                      lapply(1:length(res$m2),
                             function(x) res$m2[[x]]$ed_est)),
            type = "ed",
            dt_reps = 500,
            ci_lvl = .95,
            plot_cond = NULL,
            plot_name = NULL,
            bar_col = "#595959",
            meth_compare = meths,
            meth_sort = FALSE)
  
  ggsave(pf,
         file = "../output/graphs/exp4_ed_bias.pdf", 
         width = sp_width*2, height = sp_height*2,
         units = "cm")

# > ED CI (Facet) ####
  # Which Methods do you want to plot
  meths = rev(c("DURR_la", "IURR_la", 
                "blasso", "bridge",
                "MI_PCA",
                "MI_CART", "MI_RF", 
                "missFor", 
                "CC",
                "MI_OP",
                "GS"))
  
  pf <- plot_exp4(dt = list(lapply(1:length(res$m1),
                             function(x) res$m1[[x]]$ed_ci),
                      lapply(1:length(res$m2),
                             function(x) res$m2[[x]]$ed_ci)),
            type = "ed",
            dt_reps = 500,
            ci_lvl = .95,
            plot_cond = NULL,
            plot_name = NULL,
            bar_col = "#595959",
            meth_compare = meths,
            meth_sort = FALSE)
  
  # Save Plot
  ggsave(file = "../output/graphs/exp4_ed_ci.pdf", 
         pf,
         width = sp_width*2, height = sp_height*2,
         units = "cm")

# Time plot ---------------------------------------------------------------

# Produce data for plot
dt_outtime <- do.call(rbind, out$out_time)

# Sanitazie names
colnames(dt_outtime) <- sub("_la", "", colnames(dt_outtime))
colnames(dt_outtime) <- sub("_", "-", colnames(dt_outtime))
rownames(dt_outtime) <- c("Condition 1", "Condition 2")

# Get Plot
pf <- plot_time(
       dt = dt_outtime,
       plot_cond = NULL,
       plot_name = NULL,
       meth_compare = c(
              "DURR",
              "IURR",
              "blasso",
              "bridge",
              "MI-PCA",
              "MI-CART",
              "MI-RF",
              "stepFor",
              "MI-qp",
              "MI-am",
              "MI-OP"
       ),
       meth_sort = FALSE
)
pf

# Save Plot
ggsave(pf,
       file = "../output/graphs/exp4_time.pdf",
       width = sp_width * 2, height = sp_height * 2,
       units = "cm"
)

# Table
xtable(
       dt_outtime,
       type = "latex",
       digits = 1,
       align = c("l", rep("c", ncol(dt_outtime)))
)