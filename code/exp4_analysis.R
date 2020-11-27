### Title:    Analysis of results from experiment 2
### Project:  Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-08-27
### Notes:    reads output form results.R script and shows the numbers that
###           are used to draw the conclusions.

  rm(list = ls())
  source("./init_general.R")
  library(gridExtra)

# Read results from a run of simulation study
  # Result selection
  filename <- "exp4_simOut_20201016_2341" # old
  filename <- "exp4_simOut_20201019_1344" # current one
  filename <- "exp4_simOut_20201027_1610" # combined results from more runs
  
  # Read R object
  res <- readRDS(paste0("../output/", filename, "_res.rds"))
  out <- res
  
  # Extract names of conditions
  res$conds
  cond_names <- paste0(letters[1:nrow(res$conds)], ") ",
                       "Condition ", 1:nrow(res$conds), ": ",
                       "n = ",  res$conds$n)
  data.frame(res$conds, cond_names)
  
# Summary of set up -------------------------------------------------------
  
  # Show All
  
  # Selected paramters
  m1.indx <- c("rel")
  m2.indx <- c("NatAt")
  
  # All parameters
  m1.indx <- rownames(res$m1[[1]]$bias_per)
  m2.indx <- rownames(res$m2[[1]]$bias_per)
  
  # Bias percent
  # M1
  lapply(1:length(res$m1),
         function(x) res$m1[[x]]$bias_per[m1.indx, ])
  # M2
  lapply(1:length(res$m2),
         function(x) res$m2[[x]]$bias_per[m2.indx, ])
  
  # Bias sd
  # M1
  lapply(1:length(res$m1),
         function(x) round(res$m1[[x]]$bias_sd, 2)[m1.indx, ] )
  # M2
  lapply(1:length(res$m2),
         function(x) round(res$m2[[x]]$bias_sd, 2)[m2.indx, ] )
  
  # Confidence Intervals
  lapply(1:length(res$m1),
         function(x) res$m1[[x]]$ci_cov[m1.indx, ])[[1]]
  lapply(1:length(res$m2),
         function(x) res$m2[[x]]$ci_cov[m2.indx, ])

  
  # CI visualization
  x <- lapply(1:length(res$m1),
              function(x) res$m1[[x]]$ci_cov["rel", ])[[1]]
  x <- lapply(1:length(res$m2),
              function(x) res$m2[[x]]$ci_cov["rel", ])[[2]]
  x <- data.frame(method = factor(names(x), levels = names(x)), 
                  target = as.numeric(x-95))
  plot_target <- list(
    CI_1 = lapply(1:length(res$m1),
                  function(x) res$m1[[x]]$ci_cov["rel", ])[[1]],
    CI_2 = lapply(1:length(res$m1),
                  function(x) res$m1[[x]]$ci_cov["rel", ])[[2]]
  )
  
  lapply(plot_target, function(d){
    x <- data.frame(method = factor(names(d), levels = names(d)), 
                    target = as.numeric(d-95))
    ggplot(data = x, 
           aes(x = target, 
               y = method)) +
      # Theme First
      jtools::theme_apa() +
      # Title
      theme(plot.title = element_text(hjust = 0.5)) +
      labs(title = paste0("PROVA ", ", n = ") ,
           x     = "",
           y     = "") +
      geom_bar(stat="identity") + 
      scale_x_continuous(name = "", 
                         limits = c(-50, 5),
                         breaks = c(-50, -5, -2.5, 0, 2.5, 5),
                         labels = c(".5", ".9", ".925", ".95", ".975", "1")) +
      geom_vline(xintercept = c(-2.5, 2.5), linetype = "dashed", color = "black") +
      geom_hline(yintercept = seq(1:12), linetype = "longdash", color = "light gray") +
      coord_cartesian(xlim = c(-5, 5))
  })
  

  
# All parameters in same plot ----------------------------------------------
# This solution does not make things very clear. CI coverage is interesting,
# but not particularly interesting.
  
  bias_plots_mean <- lapply(1:2, function(p){
    plot_gg(dt = lapply(1:length(res$m2),
                        function(x) data.frame( res$m2[[x]]$bias_per))[[p]],
            parm_range = 1:13,
            type = "bias",
            y_axLab = TRUE,
            plot_name = cond_names[p],
            meth_compare = rev(c("DURR_la", "IURR_la", 
                                 "blasso", "bridge",
                                 "MI_PCA",
                                 "MI_CART", "MI_RF", 
                                 "missFor", "CC"))
    )
  } )
  
  ci_plots_mean <- lapply(1:2, function(p){
    plot_gg(dt = lapply(1:length(res$m2),
                        function(x) data.frame( res$m2[[x]]$ci_cov))[[p]],
            parm_range = 1:13,
            type = "ci",
            y_axLab = TRUE,
            plot_name = cond_names[p],
            meth_compare = rev(c("DURR_la", "IURR_la", 
                                 "blasso", "bridge",
                                 "MI_PCA",
                                 "MI_CART", "MI_RF", 
                                 "missFor", "mean"))
    )
  } )  
  
# Generic Plotting function -----------------------------------------------
# Cosutm plot_gg function for experiment 4
  plot_gg <- function(dt, 
                      dt_CIW = NULL,
                      type = "ci", 
                      plot_cond = "(empty)",
                      plot_name = "Untitled",
                      meth_compare) {
    
    # Function Inputs
    ## Example BIAS Condition 1 model 1 variable religiosity
    # dt = lapply(1:length(res$m2),
    #             function(x) res$m2[[x]]$bias_per["NatAt", ])[[2]]
    ## Example CI Condition 1, model 1, variable religiosity
    # dt = lapply(1:length(res$m2),
    #            function(x) res$m2[[x]]$ci_cov["NatAt", ])[[2]]
    # dt_CIW = lapply(1:length(res$m2),
    #             function(x) res$m2[[x]]$CIW["NatAt", ])[[2]]
    ## Example ED Condition 1, model 1 (estimates or ci?)
    # dt = lapply(1:length(res$m2),
    #            function(x) res$m2[[x]]$ed_est)[[1]]
    # dt = lapply(1:length(res$m2),
    #            function(x) res$m2[[x]]$ed_ci)[[1]]
    ## Example Confidence Interval Width
    # dt = lapply(1:length(res$m2),
    #             function(x) res$m2[[x]]$CIW["NatAt", ])[[2]]
    
    ## Generic inputs
    # type = c("bias", "ci", "ciw", "ed")[3]
    # plot_name = "Untitled"
    # plot_cond = "(empty)"
    # meth_compare = c("DURR_la", "IURR_la", "blasso", "bridge",
    #                  "MI_PCA",
    #                  "MI_CART", "MI_RF", "missFor", "CC")
    
    # Process inputs
    # var_name <- rownames(dt)
    
    if(type == "ci"){
      # Transform for graph
      # Tranform to a difference value
      x_rev <- 95 - dt
      
      # Keep only things you actualy want to show
      x_rev <- x_rev[, meth_compare]
      
      # Order by size
      rank <- order(abs(x_rev), decreasing = TRUE)
      x_rev <- x_rev[rank]
      x_rev[x_rev == 0] <- .025 
        # for display purpuses, if there are 0, 
        # substitute with small value
      
      # Make names prettier
      names(x_rev) <- sub("_la", "", names(x_rev))
      names(x_rev) <- sub("_", "-", names(x_rev))
      
      # Create CIW if requested
      if(!is.null(dt_CIW) == TRUE){
        # Keep only things you actually want to show
        CIW_plot <- as.numeric(dt_CIW[, meth_compare])
        
        # Order by size
        CIW_plot <- CIW_plot[rank]
        CIW_plot[CIW_plot == 0] <- .025
        
        # Make Character vector
        CIW_char <- paste0("(" , round(CIW_plot, 2), ")")
      } else {
        CIW_char <- ""
      }
      
      # Factorize for plot
      ggplot_input <- data.frame(method = factor(names(x_rev), 
                                                 levels = names(x_rev)), 
                                 target = as.numeric(x_rev),
                                 bar_label = CIW_char)
      
      # Limits
      plot_limits <- c(-5, 5)
      plot_breaks <- c(-5, -2.5, 0, 2.5, 5)
      plot_labels <- rev(c(".9", ".925", ".95", ".975", "1"))
      plot_vlines <- c(-2.5, 2.5)
      plot_hlines <- seq(1:12)
      plot_xlim   <- c(-5, 5)
    }

    if(type == "bias"){
      # Get rid of reference value
      x_rev <- dt[, -1]
      
      # Keep what you want to show
      x_rev <- x_rev[, meth_compare]
      
      # Order by size
      x_rev <- x_rev[order(abs(x_rev), decreasing = TRUE)]
      
      # Make names prettier
      names(x_rev) <- sub("_la", "", names(x_rev))
      names(x_rev) <- sub("_", "-", names(x_rev))
      
      # Factorize for plot
      ggplot_input <- data.frame(method = factor(names(x_rev), 
                                                 levels = names(x_rev)), 
                                 target = as.numeric(x_rev),
                                 bar_label = "")

      # Limits
      plot_xlim   <- c(-.25, .25) * 100
      plot_breaks <- c(-.25, -.1, 0, .1, .25)  * 100
      plot_labels <- c("-25", "-10", "0", "10", "25")
      plot_vlines <- c(-.1, .1)  * 100
      plot_hlines <- seq(1:nrow(ggplot_input))
    }
    
    if(type == "ciw"){
      # Keep what you want to show
      ref <- dt$GS
      x_rev <- dt[, meth_compare]
      
      # Order by size
      x_rev <- x_rev[order(abs(x_rev), decreasing = TRUE)]
      
      # Make names prettier
      names(x_rev) <- sub("_la", "", names(x_rev))
      names(x_rev) <- sub("_", "-", names(x_rev))
      
      # Factorize for plot
      ggplot_input <- data.frame(method = factor(names(x_rev), 
                                                 levels = names(x_rev)), 
                                 target = as.numeric(x_rev),
                                 bar_label = "")
      
      # Limits
      ref_max <- ggplot_input[ggplot_input$method == "CC", "target"]
      ref_min <- ggplot_input[ggplot_input$method == "missFor", "target"]
      plot_xlim   <- c(0, (ref_max*1.5))
      plot_breaks <- c(0, ref_min, ref_max, ref_max*1.5)
      plot_labels <- as.character(round(plot_breaks, 2))
      plot_vlines <- NULL
      plot_hlines <- seq(1:nrow(ggplot_input))
    }

    if(type == "ed"){
      # Keep what you want to show
      x_rev <- dt[, meth_compare]
      
      # Transform for graph
      x_rev <- sort(x_rev, decreasing = TRUE)
      
      # Make names prettier
      names(x_rev) <- sub("_la", "", names(x_rev))
      names(x_rev) <- sub("_", "-", names(x_rev))
      
      # Factorize for plot
      ggplot_input <- data.frame(method = factor(names(x_rev), 
                                                 levels = names(x_rev)), 
                                 target = as.numeric(x_rev - 0),
                                 bar_label = "")
      
      # Fixed Version
      plot_xlim   <- c(0, .5)
      plot_breaks <- seq(0, .5, 
                         length.out = 5)
      
      plot_labels <- as.character(plot_breaks)
      plot_vlines <- NULL
      plot_hlines <- NULL
    }
    
    # Produce plot
    p <- ggplot(data = ggplot_input, 
                aes(x = target, 
                    y = method)) +
      # Theme (goes first)
      jtools::theme_apa() +
      
      # Title and axis labels
      labs(title = plot_name,
           x     = "", 
           y     = "") +
      theme(plot.title = element_text(size = 12, 
                                      face = "plain", 
                                      hjust = 0.5),
            axis.text = element_text(size = 10)) +
      
      # Content
      geom_bar(stat = "identity") + 
      
      # Tweaks
      geom_text(data = ggplot_input, 
                aes(x = target + .3, 
                    y = method, 
                    label = bar_label)) + 
      scale_x_continuous(breaks = plot_breaks,
                         labels = plot_labels) +
      geom_vline(xintercept = plot_vlines, 
                 linetype = "dashed", color = "black") +
      geom_hline(yintercept = plot_hlines, 
                 size = .25,
                 color = "gray") +
      coord_cartesian(xlim = plot_xlim)
    p
    
    return(p)
  }
  
# Example use of function
  res$m1[[1]]$bias_per
  res$m2[[1]]$bias_per
  
  # Bias
  plot_gg(lapply(1:length(res$m2),
                 function(x) res$m2[[x]]$bias_per["rel", ])[[2]],
          type = "bias"
  )
  
  # Confidence Interval
  plot_gg(lapply(1:length(res$m1),
                 function(x) res$m1[[x]]$ci_cov["rel", ])[[1]],
          )
  
  # Multivariate distance bias
  plot_gg(lapply(1:length(res$m2),
                 function(x) res$m2[[x]]$ed_est)[[1]],
          type = "ed"
  )
  plot_gg(lapply(1:length(res$m2),
                 function(x) res$m2[[x]]$ed_ci)[[1]],
          type = "ed"
  )
  
# Variable of interest ----------------------------------------------------
  # Which Methods do you want to plot
  meths <- c("DURR_la", "IURR_la", 
             "blasso", "bridge",
             "MI_PCA",
             "MI_CART", "MI_RF", 
             "missFor", 
             "mean", "CC")
  
# > Bias ####
  rownames(res$m1[[1]]$bias_raw)
  rownames(res$m2[[1]]$bias_raw)

# Obtain plots per different model
  bias_plots_m1 <- lapply(1:2, function(p){
    plot_gg(lapply(1:length(res$m1),
                   function(x) res$m1[[x]]$bias_per["rel", ])[[p]], 
            type = "bias",
            plot_name = cond_names[[p]],
            meth_compare = meths) } )
  bias_plots_m2 <- lapply(1:2, function(p){
    plot_gg(lapply(1:length(res$m2),
                   function(x) res$m2[[x]]$bias_per["NatAt", ])[[p]], 
            type = "bias",
            plot_name = "",
            meth_compare = meths) } )
  
  # Display plots
  p1 <- arrangeGrob(grobs = bias_plots_m1,
                    top = "Model 1: \u03B2 Religiosity",
                    ncol = 2)
  p2 <- arrangeGrob(grobs = bias_plots_m2,
                    # top = "Percentage Relative Bias (PRB) beta Nativist Attitudes",
                    top = "Model 2: \u03B2 Nativist Attitudes",
                    ncol = 2)
  pf <- grid.arrange(p1, p2)
  
  # Save Plot
  ggsave(pf,
         file = "../output/graphs/exp4_imp_bias.pdf", 
         width = 8.25, height = 5.25,
         device = cairo_pdf)
  
# > CI   ####

  # Define plots
  ci_plots_m1 <- lapply(1:2, function(p){
    plot_gg(dt        = lapply(1:length(res$m1),
                               function(x) res$m1[[x]]$ci_cov["rel", ])[[p]],
            type      = "ci",
            plot_name = cond_names[[p]],
            meth_compare = meths) } )
  
  ci_plots_m2 <- lapply(1:2, function(p){
    plot_gg(dt        = lapply(1:length(res$m2),
                               function(x) res$m2[[x]]$ci_cov["NatAt", ])[[p]],
            type      = "ci",
            plot_name = cond_names[[p]],
            meth_compare = meths) } )
  
  # Combine and arrange plots
  p1 <- arrangeGrob(grobs = ci_plots_m1,
                    top = "Model 1: \u03B2 Religiosity", 
                    ncol = 2)
  p2 <- arrangeGrob(grobs = ci_plots_m2,
                    top = "Model 2: \u03B2 Nativist Attitudes", 
                    ncol = 2)
  pf <- grid.arrange(p1, p2)
  
   # Save Plot
  ggsave(pf, 
         file = "../output/graphs/exp4_imp_ci.pdf", 
         width = 8.25, height = 5.25,
         device = cairo_pdf)
  
# > CIW   ####
  
  # Define plots
  ci_plots_m1 <- lapply(1:2, function(p){
    plot_gg(dt        = lapply(1:length(res$m1),
                               function(x) res$m1[[x]]$CIW["rel", ])[[p]],
            type      = "ciw",
            plot_name = cond_names[[p]],
            meth_compare = meths) } )
  
  ci_plots_m2 <- lapply(1:2, function(p){
    plot_gg(dt        = lapply(1:length(res$m2),
                               function(x) res$m2[[x]]$CIW["NatAt", ])[[p]],
            type      = "ciw",
            plot_name = cond_names[[p]],
            meth_compare = meths) } )
  
  # Combine and arrange plots
  p1 <- arrangeGrob(grobs = ci_plots_m1,
                    top = "Model 1: \u03B2 Religiosity", 
                    ncol = 2)
  p2 <- arrangeGrob(grobs = ci_plots_m2,
                    top = "Model 2: \u03B2 Nativist Attitudes", 
                    ncol = 2)
  pf <- grid.arrange(p1, p2)
  
  # Save Plot
  ggsave(pf, 
         file = "../output/graphs/exp4_imp_ciw.pdf", 
         width = 8.25, height = 5.25,
         device = cairo_pdf)
  
# Multivariate distance ---------------------------------------------------
  meths <- c("DURR_la", "IURR_la", 
             "blasso", "bridge",
             "MI_PCA",
             "MI_CART", "MI_RF", 
             "missFor", 
             "mean", "CC")
# > Bias ####
  # Define plots
  bias_plots_m1 <- lapply(1:2, function(p){
    plot_gg(dt        = lapply(1:length(res$m1),
                               function(x) res$m1[[x]]$ed_est)[[p]], 
            type      = "ed",
            plot_name = cond_names[[p]],
            meth_compare = meths ) } )
  
  bias_plots_m2 <- lapply(1:2, function(p){
    plot_gg(dt        = lapply(1:length(res$m2),
                               function(x) res$m2[[x]]$ed_est)[[p]],
            type      = "ed",
            plot_name = cond_names[[p]],
            meth_compare = meths ) } )
  
  # Combine and arrange plots
  p1 <- arrangeGrob(grobs = bias_plots_m1,
                    top = "Model 1", ncol = 2
  )
  p2 <- arrangeGrob(grobs = bias_plots_m2,
                    top = "Model 2",
                    ncol = 2
  )
  pf <- grid.arrange(p1, p2)
  grid.arrange(p2) # single plot
  
  # Save Plot
  ggsave(pf,
         file = "../output/graphs/exp4_ed_bias.pdf", 
         width = 8.25, height = 5.25)
  
# > CI ####
  # Define plots
  ci_plots_m1 <- lapply(1:2, function(p){
    plot_gg(dt        = lapply(1:length(res$m1),
                               function(x) res$m1[[x]]$ed_ci)[[p]], 
            type      = "ed",
            plot_name = cond_names[[p]],
            meth_compare = meths
    ) } )
  ci_plots_m2 <- lapply(1:2, function(p){
    plot_gg(dt        = lapply(1:length(res$m2),
                               function(x) res$m2[[x]]$ed_ci)[[p]],
            type      = "ed",
            plot_name = cond_names[[p]],
            meth_compare = meths
    ) } )
  
  # Combine and arrange plots
  p1 <- arrangeGrob(grobs = ci_plots_m1,
                    top = "Model 1", ncol = 2
  )
  p2 <- arrangeGrob(grobs = ci_plots_m2,
                    top = "Model 2", ncol = 2
  )
  pf <- grid.arrange(p1, p2)
  
  # Save Plot
  ggsave(file = "../output/graphs/exp4_ed_ci.pdf", pf,
         width = 8.25, height = 5.25)
  
# > Confidence Interval Width ####
  
  res$m1[[1]]$CIW
  
  ci_res_list <- list(lapply(1:length(res$m1),
                             function(x) res$m1[[x]]$CIW),
                      lapply(1:length(res$m2),
                             function(x) res$m2[[x]]$CIW))
  lapply(1:length(ci_res_list), function(i){
    t(sapply(1:length(ci_res_list[[i]]), function(j){
      sapply(1:ncol(ci_res_list[[i]][[j]]), function(k){
        round(dist(rbind(ci_res_list[[i]][[j]][, k],
                              ci_res_list[[i]][[j]][, "GS"]),
                   method = "euclidean"), 
              3)
      })
    }))
  })
  colnames(ci_res_list[[1]][[1]])

# Time plot ---------------------------------------------------------------

  rm(list = ls())
  source("./init_general.R")
  
  # Result selection
  filename <- "exp4_simOut_20201019_1344" # current one
  filename <- "exp4_simOut_20201027_1610" # current one
  
  # Read R object
  res <- readRDS(paste0("../output/", filename, ".rds"))
  # Extract names of conditions
  res$conds
  cond_names <- paste0(letters[1:nrow(res$conds)], ") ",
                       "Condition ", 1:nrow(res$conds), ": ",
                       "n = ",  res$conds$n)
  data.frame(res$conds, cond_names)
  
  # Produce data for plot
  out_time <- sapply(1:length(names(res[[1]])), res_sem_time, out = res)
  colnames(out_time) <- names(res[[1]])
  t(out_time)
  
  # Function
  plot_gg_time <- function(dt,
                      plot_cond = "(empty)",
                      plot_name = "Untitled",
                      meth_compare) {
  ## Example Input
  # dt <- round(out_time[, 2], 1)
  # plot_cond = "(empty)"
  # plot_name = "Untitled"
  # meth_compare = c("DURR_la", "IURR_la", "blasso", "bridge",
  #                  "MI_PCA",
  #                  "MI_CART", "MI_RF")
    
  ## Prepare data for plot
  x_rev <- dt[meth_compare]# Keep what you want to show
  x_rev <- sort(x_rev, decreasing = TRUE)
  names(x_rev) <- sub("_la", "", names(x_rev))
  names(x_rev) <- sub("_", "-", names(x_rev))
  ggplot_input <- data.frame(method = factor(names(x_rev), 
                                             levels = names(x_rev)), 
                             target = as.numeric(x_rev - 0))
  
  ## Prepare plot for data
  plot_limits <- c(0, 90) # 70 min
  plot_breaks <- c(0, 30, 60, 90)
  plot_labels <- c("", "30min", "1h", "1h 30min")
  plot_vlines <- c(30, 60)
  plot_hlines <- NULL
  plot_xlim   <- c(0, 90)
  
  ## Obtain plot
  p <- ggplot(data = ggplot_input, 
              aes(x = target, 
                  y = method)) +
    # Theme (goes first)
    jtools::theme_apa() +
    
    # Title and axis labels
    labs(title = plot_name,
         x     = "", 
         y     = "") +
    theme(plot.title = element_text(size = 12, 
                                    face = "plain", 
                                    hjust = 0.5),
          axis.text = element_text(size = 10)) +
    
    # Content
    geom_bar(position = 'dodge', stat = "identity") + 
    geom_text(aes(label = target), hjust = -.2) + 
    
    # Tweaks
    scale_x_continuous(breaks = plot_breaks,
                       labels = plot_labels) +
    geom_vline(xintercept = plot_vlines, 
               linetype = "dashed", color = "black") +
    geom_hline(yintercept = plot_hlines, 
               size = .25,
               color = "gray") +
    coord_cartesian(xlim = plot_xlim)
  
  return(p)
  
  }
  
  ## Use fucntion on time data
  p <- lapply(1:ncol(out_time), function(j){
    plot_gg_time( round(out_time[, j],1),
                  plot_name = cond_names[j],
                  meth_compare = c("DURR_la", "IURR_la", "blasso", "bridge",
                                   "MI_PCA",
                                   "MI_CART", "MI_RF")
                  )
  })
  # Reuglar Presentation
  p_time <- arrangeGrob(grobs = p,
                   top = "Average Imputation time (minutes) for each MI method",
                   ncol = 2)
  # Poster presentation
  p_time <- arrangeGrob(grobs = p,
                        ncol = 2)
  # Plot
  pf <- grid.arrange(p_time)

  ggsave(file = "../output/graphs/exp4_time.pdf", pf,
         width = 8.25, height = 5.25/2)
