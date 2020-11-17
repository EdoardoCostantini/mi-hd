### Title:    Analysis of results
### Project:  Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-07-09
### Notes:    reads output form results.R script and shows the numbers that
###           are used to draw the conclusions.

  rm(list = ls())
  source("./init_general.R")

# Read results from a run of simulation study
  # Checks
  res <- exp1_res <- readRDS("../output/exp1_simOut_20200731_1735_res.rds")
  res <- exp1_res <- readRDS("../output/exp1_simOut_20200801_1620_res.rds")
  # They are equivalent, but second file has more repetitions (500 vs 750)
  
# Bias --------------------------------------------------------------------

  # Result 1 - Easy round, all good except CART and RF
  sum_exp1_sem$cond1$cond
  sum_exp1_sem$cond1$bias_per[c(1,7,13), ]
  
  # Result 2 - Covariance trouble for Blasso and DURR w/ p 50 pm .3
  sum_exp1_sem$cond3$cond
  sum_exp1_sem$cond3$bias_per[c(1,7,13), ]
  sinfo <- session_info()
  ls(sinfo)
  sinfo$platform
  list(Platform = sinfo$platform,
       Pacakges = data.frame(sinfo$packages))
  # Result 3 - Variance and covariance trouble for MI-PCA and IURR, respectively
  sum_exp1_sem$cond4$cond
  sum_exp1_sem$cond4$bias_per[c(1,7,13), ]

# CI ----------------------------------------------------------------------
  
  # Result 2 - Coverage is a function of pm, not p (report mean CI coverage)
  # pm .1
  t(round(sapply(list(p50 = sum_exp1_sem$cond1$ci_cov,  # p 50 
                      p500= sum_exp1_sem$cond2$ci_cov), # p 500
               colMeans), 0))
  # p 50
  t(round(sapply(list(pm.1 = sum_exp1_sem$cond1$ci_cov,  # pm .1
                      pm.3 = sum_exp1_sem$cond3$ci_cov), # pm .3
                 colMeans), 0))
  
  # p 500
  t(round(sapply(list(pm.1 = sum_exp1_sem$cond2$ci_cov,  # pm .1
                      pm.3 = sum_exp1_sem$cond4$ci_cov), # pm .3
                 colMeans), 0))
  
# Multivairate Distance Analysis ------------------------------------------

  # Sem model estimates 
  # All parameters
  sem_ed_all <- round(t(sapply(list(p50pm.1  = sum_exp1_sem$cond1, 
                                    p500pm.1 = sum_exp1_sem$cond2, 
                                    p50pm.3  = sum_exp1_sem$cond3, 
                                    p500pm.3 = sum_exp1_sem$cond4), 
                               res_ed_est, 
                               measure = "all"
  )), 3)
  
  # Means
  sem_ed_mu <- round(t(sapply(list(p50pm.1  = sum_exp1_sem$cond1, 
                                   p500pm.1 = sum_exp1_sem$cond2, 
                                   p50pm.3  = sum_exp1_sem$cond3, 
                                   p500pm.3 = sum_exp1_sem$cond4), 
                              res_ed_est, 
                              measure = "mean"
  )), 3)
  
  # Variances
  sem_ed_var <- round(t(sapply(list(p50pm.1  = sum_exp1_sem$cond1, 
                                    p500pm.1 = sum_exp1_sem$cond2, 
                                    p50pm.3  = sum_exp1_sem$cond3, 
                                    p500pm.3 = sum_exp1_sem$cond4), 
                               res_ed_est, 
                               measure = "var"
  )), 3)
  
  # Covariances
  sem_ed_cov <- round(t(sapply(list(p50pm.1  = sum_exp1_sem$cond1, 
                                    p500pm.1 = sum_exp1_sem$cond2, 
                                    p50pm.3  = sum_exp1_sem$cond3, 
                                    p500pm.3 = sum_exp1_sem$cond4), 
                               res_ed_est, 
                               measure = "cov"
  )), 3)
  round(rbind(sem_ed_all, 
              sem_ed_mu, 
              sem_ed_var, 
              sem_ed_cov), 2)
  
  round(sem_ed_mu, 2)
  round(sem_ed_var, 2)
  round(sem_ed_cov, 2)
  
  # Make latex tables
  latex_tab_input <- rbind(sem_ed_all, rep(NA, ncol(sem_ed_all)),
                           sem_ed_mu, rep(NA, ncol(sem_ed_all)),
                           sem_ed_var, rep(NA, ncol(sem_ed_all)),
                           sem_ed_cov)
  latex_tab_input <- do.call(rbind,
                             lapply(list(sem_ed_all, 
                                         sem_ed_mu, 
                                         sem_ed_var, 
                                         sem_ed_cov), 
                                    function(x){
                                      rownames(x) <- c("p = 50, pm = .1",
                                                       "p = 50, pm = .3",
                                                       "p = 500, pm = .1",
                                                       "p = 500, pm = .3")
                                      x}
                             )
  )
  xtable(latex_tab_input,
         align = c("l", rep("c", ncol(latex_tab_input))))
  
  # or
  lapply(list(sem_ed_all, sem_ed_mu, sem_ed_var, sem_ed_cov), function(x){
    rownames(x) <- c("p = 50, pm = .1",
                     "p = 50, pm = .3",
                     "p = 500, pm = .1",
                     "p = 500, pm = .3")
    xtable(x,
           align = c("l", rep("c", ncol(x))) )
  }
  )
  

# Summary Table -----------------------------------------------------------
# This is a selected paramters version for summary paper presentation.

  col_id <- colnames(sum_exp1_sem$cond4$bias_raw)[c(12, 10, 2:9)]
  indx <- rownames(sum_exp1_sem$cond4$bias_raw)[c(1, 4, 
                                                  7, 10, 
                                                  13, 15, 25)]
  sum_exp1_sem$cond4$cond
  sum_exp1_sem$cond4$ci_cov[indx, ]
  
  sum_exp1_sem$cond4$cond
  sum_exp1_sem$cond4$bias_per[indx, ]
  
  genTableEAM <- function(x){
    # Generates section of the table for EAM turn in paper
    store <- NULL
    for (i in 1:length(indx)) {
      store <- cbind(store,
                     t(sum_exp1_sem[[x]]$bias_per[indx, col_id])[, i],
                     t(sum_exp1_sem[[x]]$ci_cov[indx, col_id])[, i])
    }
    return(store)
  }
  
  list_tables <- lapply(list(cond1=1,
                             cond2=2,
                             cond3=3,
                             cond4=4), 
                        genTableEAM)
  
  mt_table <- do.call(rbind, lapply(list_tables, rbind, NA))

  write.csv(mt_table, "~/Desktop/mt_table.csv")
  
  xtable(mt_table,
         align = c("l", rep("c", ncol(mt_table))))
  
# LM models ---------------------------------------------------------------

  # Bias
  sum_exp1_lm$cond1$cond
  sum_exp1_lm$cond1$bias_per
  sum_exp1_lm$cond3$cond
  sum_exp1_lm$cond3$bias_per
  
    # HighD
  sum_exp1_lm$cond2$cond
  sum_exp1_lm$cond2$bias_per
  sum_exp1_lm$cond4$cond
  sum_exp1_lm$cond4$bias_per
  
  # CI
  sum_exp1_lm$cond1$cond
  sum_exp1_lm$cond1$ci_cov
  sum_exp1_lm$cond3$cond
  sum_exp1_lm$cond3$ci_cov
  
    # HighD
  sum_exp1_lm$cond2$cond
  sum_exp1_lm$cond2$ci_cov
  sum_exp1_lm$cond4$cond
  sum_exp1_lm$cond4$ci_cov
  

# Plots -------------------------------------------------------------------

  plot_gg <- function(dt,
                      type = "bias", 
                      plot_cond = "(empty)",
                      plot_name = NULL,
                      parm_range = 1:2,
                      y_axLab = TRUE,
                      meth_compare) {
    ## Function inputs
    ## Generic
    # y_axLab = TRUE # say I want the labes
    # parm_range = 1:6
    # type = "ci"
    # plot_name = "Untitled"
    # meth_compare = rev(c("DURR_la", "IURR_la", "blasso", "bridge",
    #                      "MI_PCA",
    #                      "MI_CART", "MI_RF", "missFor", "CC"))
    ## Bias
    # dt = lapply(1:length(res$sem),
    #             function(x) data.frame( res$sem[[x]]$bias_per))[[4]]
    ## CIR
    # dt = lapply(1:length(res$sem),
    #            function(x) res$sem[[x]]$ci_cov)[[1]]
    
    ## Prep data for plot
    # Select range of interest
    dt_edit <- dt[parm_range, ]
    
    # Order Methods
    dt_edit <- dt_edit[, meth_compare]
    
    # Make names more pretty
    colnames(dt_edit) <- sub("_la", "", colnames(dt_edit))
    colnames(dt_edit) <- sub("_", "-", colnames(dt_edit))
    
    # Shape for ggplot
    dt_edit[nrow(dt_edit)+1,] <- 0  # add blank row to improve visualization
    dt_edit <- dt_edit %>% gather()
    dt_edit$id <- 1:nrow(dt_edit)  # add an 
    
    # Ticks 
    if(type == "bias"){
      plot_xlim   <- c(-20, 20)
      plot_xbreaks <- c(-20, -10, 0, 10, 20)
      plot_xlabels <- rev(c("-20", "-10", "0", "10", "20"))
      
      step_size   <- max(dt_edit$id)/length(unique(dt_edit$key)) / 2
      plot_steps  <- seq(0, nrow(dt_edit), by = step_size)
      plot_ybreaks <- plot_steps[c(FALSE, TRUE)] # keep every other element
      if(y_axLab == TRUE){
        plot_ylabels <- as.character(unique(dt_edit$key))
      } else {
        plot_ylabels <- NULL
      }
      plot_vlines <- plot_steps[c(TRUE, FALSE)] # keep every other element
      plot_hlines <- c(-10, 10)
    }
    
    if(type == "bias_raw"){
      plot_xlim   <- c(-1, 1)
      plot_xbreaks <- c(-1, -.5, 0, .5, 1)
      plot_xlabels <- rev(c("-1", "-.5", "0", ".5", "1"))
      
      step_size    <- max(dt_edit$id)/length(unique(dt_edit$key)) / 2
      plot_steps   <- seq(0, nrow(dt_edit), by = step_size)
      plot_ybreaks <- plot_steps[c(FALSE, TRUE)] # keep every other element
      if(y_axLab == TRUE){
        plot_ylabels <- as.character(unique(dt_edit$key))
      } else {
        plot_ylabels <- NULL
      }
      plot_vlines <- plot_steps[c(TRUE, FALSE)] # keep every other element
      plot_hlines <- c(-.5, .5)
    }
    
    if(type == "ci"){
      dt_edit$value[c(rep(TRUE, length(parm_range)), FALSE)] <- 
        95 - dt_edit$value[c(rep(TRUE, length(parm_range)), FALSE)]
      
      plot_limits <- c(-5, 15)
      
      plot_xlim   <- c(-5, 15)
      plot_xbreaks <- c(-5, -2.5, 0, 2.5, 5, 15)
      plot_xlabels <- rev(c(".8", ".9", ".925", ".95", ".975", "1"))
      
      step_size   <- max(dt_edit$id)/length(unique(dt_edit$key)) / 2
      plot_steps  <- seq(0, nrow(dt_edit), by = step_size)
      plot_ybreaks <- plot_steps[c(FALSE, TRUE)] # keep every other element
      if(y_axLab == TRUE){
        plot_ylabels <- as.character(unique(dt_edit$key))
      } else {
        plot_ylabels <- NULL
      }
      plot_vlines <- plot_steps[c(TRUE, FALSE)] # keep every other element
      plot_hlines <- c(-2.5, 2.5)
    }
    
    # Plot
    p <- ggplot(dt_edit, aes(x = value, y = id)) +
      # Theme (goes first)
      jtools::theme_apa() +
      
      # Title and axis labels
      labs(title = plot_name,
           x     = "", 
           y     = "") +
      theme(plot.title = element_text(hjust = 0.5),
            axis.title = element_text(size = 5),
            axis.text = element_text(size = 5)) +
      
      # Content
      geom_segment(aes(xend = 0, 
                       yend = id),
                   color = "gray") + 
      
      # Tweaks
      scale_y_continuous(breaks = plot_ybreaks,
                         labels = plot_ylabels) +
      scale_x_continuous(breaks = plot_xbreaks,
                         labels = plot_xlabels) +
      geom_vline(xintercept = plot_hlines,
                 linetype = "dashed", color = "black") +
      geom_hline(yintercept = plot_vlines,
                 size = .25, 
                 color = "black") +
      coord_cartesian(xlim = plot_xlim)
    p
    return(p)
  }

  ## > Bias (Single Plot) #####
  plot_gg(dt = lapply(1:length(res$sem),
                      function(x) data.frame( res$sem[[x]]$bias_per))[[4]],
          parm_range = 1:6,
          type = "bias",
          meth_compare = rev(c("DURR_la", "IURR_la", 
                               "blasso", "bridge",
                               "MI_PCA",
                               "MI_CART", "MI_RF", 
                               "missFor", "CC"))
  )
  
  ## > Bias (Grid distinguishing condition and measure) #####
  bias_plots_mean <- lapply(1:4, function(p){
    plot_gg(dt = lapply(1:length(res$sem),
                        function(x) data.frame( res$sem[[x]]$bias_per))[[p]],
            parm_range = 1:6,
            type = "bias",
            meth_compare = rev(c("DURR_la", "IURR_la", 
                                 "blasso", "bridge",
                                 "MI_PCA",
                                 "MI_CART", "MI_RF", 
                                 "missFor", "CC"))
    )
  } )  
  bias_plots_var <- lapply(1:4, function(p){
    plot_gg(dt = lapply(1:length(res$sem),
                        function(x) data.frame( res$sem[[x]]$bias_per))[[p]],
            parm_range = 7:12,
            type = "bias",
            meth_compare = rev(c("DURR_la", "IURR_la", 
                                 "blasso", "bridge",
                                 "MI_PCA",
                                 "MI_CART", "MI_RF", 
                                 "missFor", "CC"))
    )
  } )  
  bias_plots_cov <- lapply(1:4, function(p){
    plot_gg(dt = lapply(1:length(res$sem),
                        function(x) data.frame( res$sem[[x]]$bias_per))[[p]],
            parm_range = -(1:12),
            type = "bias",
            meth_compare = rev(c("DURR_la", "IURR_la", 
                                 "blasso", "bridge",
                                 "MI_PCA",
                                 "MI_CART", "MI_RF", 
                                 "missFor", "CC"))
    )
  } )  
  
  # Display Plots
  p1 <- arrangeGrob(grobs = bias_plots_mean,
                    top = "Mean", 
                    ncol = 1)
  p2 <- arrangeGrob(grobs = bias_plots_var,
                    top = "Variances",
                    ncol = 1)
  p3 <- arrangeGrob(grobs = bias_plots_cov,
                    top = "Covariances",
                    ncol = 1)
  pf <- grid.arrange(p1, p2, p3, ncol = 3)
  pf
  ggsave(file  = "~/Desktop/bias.pdf", 
         width = 8.25, height = 11.75,
         arrangeGrob(p1, p2, p3, ncol = 3))

  ## > CI (Single Plot) #####
  plot_gg(dt = lapply(1:length(res$sem),
                      function(x) data.frame( res$sem[[x]]$ci_cov))[[1]],
          parm_range = 1:6,
          type = "ci",
          meth_compare = rev(c("DURR_la", "IURR_la", 
                               "blasso", "bridge",
                               "MI_PCA",
                               "MI_CART", "MI_RF", 
                               "missFor", "CC"))
  )
  
  ## > CI (Grid distinguishing condition and measure) #####
  ci_plots_mean <- lapply(1:4, function(p){
    plot_gg(dt = lapply(1:length(res$sem),
                        function(x) data.frame( res$sem[[x]]$ci_cov))[[p]],
            parm_range = 1:6,
            type = "ci",
            meth_compare = rev(c("DURR_la", "IURR_la", 
                                 "blasso", "bridge",
                                 "MI_PCA",
                                 "MI_CART", "MI_RF", 
                                 "missFor", "CC"))
    )
  } )  
  ci_plots_var <- lapply(1:4, function(p){
    plot_gg(dt = lapply(1:length(res$sem),
                        function(x) data.frame( res$sem[[x]]$ci_cov))[[p]],
            parm_range = 7:12,
            type = "ci",
            y_axLab = TRUE,
            meth_compare = rev(c("DURR_la", "IURR_la", 
                                 "blasso", "bridge",
                                 "MI_PCA",
                                 "MI_CART", "MI_RF", 
                                 "missFor", "CC"))
    )
  } )  
  ci_plots_cov <- lapply(1:4, function(p){
    plot_gg(dt = lapply(1:length(res$sem),
                        function(x) data.frame( res$sem[[x]]$ci_cov))[[p]],
            parm_range = -(1:12),
            type = "ci",
            y_axLab = TRUE,
            meth_compare = rev(c("DURR_la", "IURR_la", 
                                 "blasso", "bridge",
                                 "MI_PCA",
                                 "MI_CART", "MI_RF", 
                                 "missFor", "CC"))
    )
  } )  
  
  # Display Plots
  p1 <- arrangeGrob(grobs = ci_plots_mean,
                    top = "Mean", 
                    ncol = 1)
  p2 <- arrangeGrob(grobs = ci_plots_var,
                    top = "Variances",
                    ncol = 1)
  p3 <- arrangeGrob(grobs = ci_plots_cov,
                    top = "Covariances",
                    ncol = 1)
  pf <- grid.arrange(p1, p2, p3, ncol = 3)
  pf
  ggsave(file  = "~/Desktop/ci.pdf", 
         width = 8.25, height = 11.75,
         arrangeGrob(p1, p2, p3, ncol = 3))