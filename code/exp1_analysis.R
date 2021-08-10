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
  # Full information names
  cond_parms <- paste0("p = ",  res$conds$p, "; ",
                       "pm = ",  res$conds$pm)
  cond_labels <- c("low-dim-low-pm",
                  "high-dim-low-pm",
                  "low-dim-high-pm",
                  "high-dim-high-pm")
  cond_names <- paste(cond_labels, cond_parms, sep = " \n ")
  
  # Numbered Conditions
  # cond_names <- paste0("Condition ", 1:4)
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
                cond_labels = cond_names,
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
                cond_labels = cond_names,
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
                cond_labels = cond_names,
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
                cond_labels = cond_names,
                ci_lvl = .95,
                meth_compare = meth_compare)
  pf
  ggsave(file = "../output/graphs/exp1_CI.pdf",
         width = gp_width, height = gp_height,
         units = "cm",
         pf)

# New version -------------------------------------------------------------

  output_sem <- gg_out_sem

  # Condition Names / Labels
  label_cond <- unique(output_sem$cond)
  label_parm <- unique(output_sem$parm)

  # Threshold lines
  vline_bias <- 10

  # SE for threshold
  ci_lvl <- .95
  dt_reps <- 1e3
  SEp <- sqrt(ci_lvl*(1-ci_lvl)/dt_reps)
  low_thr <- (.95-SEp*2)*100
  hig_thr <- (.95+SEp*2)*100
  vline_burton <- c(low_thr, hig_thr)
  vline_vanBuu <- 90

  # Bias
  x <- 1
  methods_sel <- levels(output_sem$methods)[1:9]
  pm_grep <- "0.1"

  pf <- output_sem %>%
    filter(analysis == unique(analysis)[x],
           grepl(pm_grep, cond),
           variable %in% c("Min", "Mean", "Max"),
           methods %in% methods_sel
    ) %>%
    # Drop pm = ** as we are plotting only one value for this condition
    mutate(cond = fct_relabel(cond,
                              str_replace,
                              " pm = [0-9]\\.[0-9]", "")
    ) %>%
    mutate(methods = fct_relabel(methods,
                                 str_replace,
                                 "-la", "")
    ) %>%
    # Main Plot
    ggplot(data = ., aes(y = methods,
                         x = value,
                         shape = variable)) +
    geom_point(size = 1) +
    geom_line(aes(group = methods),
              size = .25) +
    # Grid
    facet_grid(rows = vars(factor(parm,
                  levels = unique(parm))),
               cols = vars(cond)) +
    geom_vline(data = data.frame(xint = c(10),
                                 analysis = "Percentage Relative Bias"),
               linetype = "dashed",
               size = .35,
               aes(xintercept = xint)) +

    # Format
    scale_y_discrete(limits = rev) +
    scale_shape_manual(values = c(1, 4, 16)) +
    coord_cartesian(xlim = c(0, 50)) +
    labs(#title = label_parm[x],
      x     = NULL,
      y     = NULL,
      linetype = NULL,
      shape = NULL) +
    theme(legend.position ="bottom",
          text = element_text(size = 9))

    pf

    plot_name <- paste0("../output/graphs/exp1_", "bias", "_",
                        as.numeric(pm_grep)*100, ".pdf")
    ggsave(file = plot_name,
           width = gp_width, height = gp_height*.65,
           units = "cm",
           pf)

  # Confidence Interval
  x <- 2
  methods_sel <- levels(output_sem$methods)#[-8]

  # SE for threshold
  ci_lvl <- .95
  dt_reps <- 1e3
  SEp <- sqrt(ci_lvl*(1-ci_lvl)/dt_reps)
  low_thr <- (.95-SEp*2)*100
  hig_thr <- (.95+SEp*2)*100
  vline_burton <- c(low_thr, hig_thr)
  vline_vanBuu <- 90
  xci_breaks <- c(80, 90, round(vline_burton, 1), 100)

  pf <- output_sem %>%
    filter(analysis == unique(analysis)[x],
           grepl(pm_grep, cond),
           variable %in% c("Min", "Mean", "Max"),
           methods %in% methods_sel
    ) %>%
    # Drop pm = ** as we are plotting only one value for this condition
    mutate(cond = fct_relabel(cond,
                              str_replace,
                              " pm = [0-9]\\.[0-9]", "")
    ) %>%
    mutate(methods = fct_relabel(methods,
                                 str_replace,
                                 "-la", "")
    ) %>%
    # Main Plot
    ggplot(data = ., aes(y = methods,
                         x = value,
                         shape = variable)) +
    geom_point(size = 1) +
    geom_line(aes(group = methods),
              size = .25) +
    # Grid
    facet_grid(rows = vars(factor(parm,
                  levels = unique(parm))),
               cols = vars(cond)) +
    geom_vline(data = data.frame(xint = c(vline_burton),
                                     analysis = "CI coverage"),
                   aes(xintercept = xint,
                       lty = paste0("Burton's range")),
                   size = .5) +
    geom_vline(data = data.frame(xint = c(vline_vanBuu),
                                 analysis = "CI coverage"),
               linetype = "dashed",
               aes(xintercept = xint),
               size = .35) +

    # Format
    scale_y_discrete(limits = rev) +
    scale_x_continuous(labels = as.character(xci_breaks), breaks = xci_breaks) +
    coord_cartesian(xlim = c(min(xci_breaks), max(xci_breaks))) +
    scale_shape_manual(values = c(1, 4, 16)) +
    scale_linetype_manual(values = c("dotted", "dotted")) +
    labs(#title = label_parm[x],
      x     = NULL,
      y     = NULL,
      linetype = NULL,
      shape = NULL) +
    theme(legend.position ="bottom",
          text = element_text(size = 9))

    pf

    plot_name <- paste0("../output/graphs/exp1_", "CI", "_",
                        as.numeric(pm_grep)*100, ".pdf")
    ggsave(file = plot_name,
           width = gp_width, height = gp_height*.65,
           units = "cm",
           pf)

