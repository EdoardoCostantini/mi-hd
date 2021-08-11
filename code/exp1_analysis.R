### Title:    Analysis of results
### Project:  Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-07-09
### Notes:    reads output form results.R script and shows the numbers that
###           are used to draw the conclusions.

  rm(list = ls())
  source("./init_general.R")
  
# Read results from the combined results
  res <- exp1_res <- readRDS("../output/exp1_simOut_2020121516_res.rds")
  
# Plot Sizes Decisions

  gp_width <- 15
  gp_height <- 20
  sp_width <- 3.75
  sp_height <- 6

# New version -------------------------------------------------------------

  output_sem <- res

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
  pm_grep <- "0.3"

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
  vline_vanBuu <- c(90, 99)
  xci_breaks <- sort(c(80, vline_vanBuu, round(vline_burton, 1), 100))

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
               linetype = "longdash",
               aes(xintercept = xint),
               size = .35) +

    # Format
    scale_y_discrete(limits = rev) +
    scale_x_continuous(labels = as.character(xci_breaks/100),
                       breaks = xci_breaks) +
    coord_cartesian(xlim = c(min(xci_breaks), max(xci_breaks))) +
    scale_shape_manual(values = c(1, 4, 16)) +
    scale_linetype_manual(values = c("dotted", "dotted")) +
    labs(#title = label_parm[x],
      x     = NULL,
      y     = NULL,
      linetype = NULL,
      shape = NULL) +
    theme(legend.position ="bottom",
          text = element_text(size = 9),
          axis.text.x = element_text(angle = 90,
                                     vjust = 0.5,
                                     hjust = 1))

    pf

    plot_name <- paste0("../output/graphs/exp1_", "CI", "_",
                        as.numeric(pm_grep)*100, ".pdf")
    ggsave(file = plot_name,
           width = gp_width, height = gp_height*.65,
           units = "cm",
           pf)

