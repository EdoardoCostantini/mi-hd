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
  

# Generic Plotting function -----------------------------------------------
  plot_gg <- function(dt, type = "ci", 
                      plot_cond = "(empty)",
                      plot_name = "Untitled",
                      meth_compare) {
    
    # Function Inputs
    ## Example BIAS Condition 1 model 1 variable religiosity
    # dt = lapply(1:length(res$m2),
    #             function(x) res$m2[[x]]$bias_per["NatAt", ])[[2]]
    ## Example bias_worst Condition 1, model 1 (estimates or ci?)
    # dt = lapply(1:length(res$m2),
    #            function(x) res$m2[[x]]$bias_sd)[[1]]
    ## Example CI Condition 1, model 1, variable religiosity
    # dt = lapply(1:length(res$m2),
    #            function(x) res$m2[[x]]$ci_cov["NatAt", ])[[2]]
    ## Example CI worst Condition 1, model 1, variable religiosity
    # dt = lapply(1:length(res$m2),
    #            function(x) res$m2[[x]]$ci_cov["NatAt", ])[[1]]
    ## Example ED Condition 1, model 1 (estimates or ci?)
    # dt = lapply(1:length(res$m2),
    #            function(x) res$m2[[x]]$ed_est)[[1]]
    # dt = lapply(1:length(res$m2),
    #            function(x) res$m2[[x]]$ed_ci)[[1]]
    ## Example ED for condition 1 all parameters
    # dt = res$ed_all[[1]]$ed_est
    # dt = res$ed_all[[1]]$ed_ci
    
    ## Generic inputs
    # type = c("bias", "ci", "ed")[3]
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
      # extra <- extra[meth_compare]
      
      # Order by size
      # extra <- extra[order(abs(x_rev), decreasing = TRUE)]
      x_rev <- x_rev[order(abs(x_rev), decreasing = TRUE)]
      x_rev[x_rev == 0] <- .025 
        # for display purpuses, if there are 0, 
        # substitute with small value
      
      # Make names prettier
      names(x_rev) <- sub("_la", "", names(x_rev))
      names(x_rev) <- sub("_", "-", names(x_rev))
      
      # Factorize for plot
      ggplot_input <- data.frame(method = factor(names(x_rev), 
                                                 levels = names(x_rev)), 
                                 target = as.numeric(x_rev))
      
      # Limits
      plot_limits <- c(-5, 5)
      plot_breaks <- c(-5, -2.5, 0, 2.5, 5)
      plot_labels <- rev(c(".9", ".925", ".95", ".975", "1"))
      plot_vlines <- c(-2.5, 2.5)
      plot_hlines <- seq(1:12)
      plot_xlim   <- c(-5, 5)
    }
    
    if(type == "ci_worst"){
      # Transform for graph
      # Tranform to a difference value
      x <- lapply(names(dt), function(j){
        # j <- names(dt)[1]
        max.out <- dt[which.max(abs(95-dt[[j]])), ][j]
        data.frame(method  = names(max.out),
                   coef    = rownames(max.out),
                   cic     = 95 - as.numeric(max.out))
      })
      x <- do.call(rbind, x)
      
      # Keep what you want to show
      x_rev <- x[x$method %in% meth_compare, ]
      
      # Order by size
      x_rev <- x_rev[order(abs(x_rev$cic), decreasing = TRUE), ]
      
      # Make names prettier
      x_rev$method <- sub("_la", "", x_rev$method)
      x_rev$method <- sub("_", "-", x_rev$method)
      
      # Factorize for plot
      ggplot_input <- data.frame(method = factor(x_rev$method, 
                                                 levels = x_rev$method), 
                                 target = as.numeric(x_rev$cic))
      
      # ggplot_input <- ggplot_input[-nrow(ggplot_input), ]
      
      # Limits
      plot_limits <- c(-5, 50)
      plot_breaks <- c(-5, -2.5, 0, 2.5, 5, 50)
      plot_labels <- rev(c(".5", ".9", ".925", ".95", ".975", "1"))
      plot_vlines <- c(-2.5, 2.5)
      plot_hlines <- seq(1:12)
      plot_xlim   <- c(-5, 50)
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
                                 target = as.numeric(x_rev))
      
      # ggplot_input <- ggplot_input[-nrow(ggplot_input), ]
      
      # Limits
      plot_xlim   <- c(-.25, .25) * 100
      plot_breaks <- c(-.25, -.1, 0, .1, .25)  * 100
      plot_labels <- c("-25", "-10", "0", "10", "25")
      plot_vlines <- c(-.1, .1)  * 100
      plot_hlines <- seq(1:nrow(ggplot_input))
    }
    
    if(type == "bias_worst"){
      # Transform for graph
      # Tranform to a difference value
      x <- lapply(names(dt), function(j){
        # j <- names(dt)[1]
        max.out <- dt[which.max(abs(dt[[j]])), ][j]
        data.frame(method  = names(max.out),
                   coef    = rownames(max.out),
                   bias_sd = as.numeric(max.out))
      })
      x <- do.call(rbind, x)
      
      # Keep what you want to show
      x_rev <- x[x$method %in% meth_compare, ]
      
      # Order by size
      x_rev <- x_rev[order(abs(x_rev$bias_sd), decreasing = TRUE), ]
      
      # Make names prettier
      x_rev$method <- sub("_la", "", x_rev$method)
      x_rev$method <- sub("_", "-", x_rev$method)
      
      # Factorize for plot
      ggplot_input <- data.frame(method = factor(x_rev$method, 
                                                 levels = x_rev$method), 
                                 target = as.numeric(x_rev$bias_sd))
      
      # ggplot_input <- ggplot_input[-nrow(ggplot_input), ]
      
      # Limits
      plot_xlim   <- c(-2, 2)
      plot_breaks <- c(-2, -1, 0, 1, 2)
      plot_labels <- c("-2", "-1", "0", "1", "2")
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
      
      # # Separate GS from the bunch
      # x_rev_GS <- x_rev["GS"]
      # x_rev <- x_rev[-which(names(x_rev)=="GS")]
      # 
      # # Include GS as the first in comparison
      # x_rev <- cbind(x_rev[order(abs(x_rev), decreasing = TRUE)], 
      #                x_rev_GS)
      
      # Factorize for plot
      ggplot_input <- data.frame(method = factor(names(x_rev), 
                                                 levels = names(x_rev)), 
                                 target = as.numeric(x_rev - 0))
      
      # # Get rid of GS if it is a bias plot (not intersting that it is 0)
      # GS_index <- which(ggplot_input$method == "GS")
      # if(ggplot_input[GS_index, "target"] == 0){
      #   ggplot_input <- ggplot_input[-GS_index, ]
      # } 
      
      # Limits
      # Flexible version
      # plot_xlim   <- c(floor(min(ggplot_input$target)), 
      #                  ceiling(max(ggplot_input$target)))
      # plot_breaks <- seq(floor(min(ggplot_input$target)), 
      #                    ceiling(max(ggplot_input$target)), 
      #                    length.out = 5)
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
      theme(plot.title = element_text(hjust = 0.5),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 12)) +
      
      # Content
      geom_bar(stat = "identity") + 
      
      # Tweaks
      scale_x_continuous(breaks = plot_breaks,
                         labels = plot_labels) +
      geom_vline(xintercept = plot_vlines, 
                 linetype = "dashed", color = "black") +
      geom_hline(yintercept = plot_hlines, 
                 size = .25, 
                 # linetype = "longdash", 
                 color = "gray") +
      coord_cartesian(xlim = plot_xlim)
    p
    return(p)
  }
  
  
# Example use of function
  res$m1[[1]]$bias_per
  res$m2[[1]]$bias_per
  res$m2[[1]]
  
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
  
# > Bias ####
  rownames(res$m1[[1]]$bias_raw)
  rownames(res$m2[[1]]$bias_raw)
  # Extract what you want to plot
  bias_res_list <- list(lapply(1:length(res$m1),
                               function(x) res$m1[[x]]$bias_per["rel", ]),
                        lapply(1:length(res$m2),
                               function(x) res$m2[[x]]$bias_per["NatAt", ]))
  
  # Obtain plots per different condition
  bias_plots_c1 <- lapply(1:length(bias_res_list), function(p){
    plot_gg(bias_res_list[[p]][[1]], 
            type = "bias",
            plot_name = c("Low Dimensional Condition"),
            meth_compare = c("DURR_la", "IURR_la", "blasso", "bridge",
                             "MI_PCA",
                             "MI_CART", "MI_RF", "mean", "CC")
            )
  } )
  bias_plots_c2 <- lapply(1:length(bias_res_list), function(p){
    plot_gg(bias_res_list[[p]][[2]], 
            type = "bias",
            plot_name = c("High Dimensional Condition"),
            meth_compare = c("DURR_la", "IURR_la", "blasso", "bridge",
                             "MI_PCA",
                             "MI_CART", "MI_RF", "mean", "CC")
    )
  } )
  
  # Display plots
  p1 <- arrangeGrob(grobs = list(bias_plots_c1[[1]], bias_plots_c2[[1]]),
                    top = "Model 1: \u03B2 Religiosity", ncol = 2,
  )
  p2 <- arrangeGrob(grobs = list(bias_plots_c1[[2]], bias_plots_c2[[2]]),
                    top = "Percentage Relative Bias (PRB) Nativist Attitudes' \u03B2",
                    ncol = 2
  )
  pf <- grid.arrange(p1, p2)
  pf <- grid.arrange(p2) # single plot
  
  ggsave(file = "../output/graphs/exp4_imp_bias.pdf", pf, device = cairo_pdf)
  
  ## Wrost Bias
  
  # Extract what you want to plot
  bias_res_list <- list(lapply(1:length(res$m1),
                               function(x) res$m1[[x]]$bias_sd),
                        lapply(1:length(res$m2),
                               function(x) res$m2[[x]]$bias_sd))
  
  # Obtain plots per different condition
  bias_plots_c1 <- lapply(1:length(bias_res_list), function(p){
    plot_gg(bias_res_list[[p]][[1]], 
            type = "bias_worst",
            plot_name = c("Low Dimensional Condition"),
            meth_compare = c("DURR_la", "IURR_la", "blasso", "bridge",
                             "MI_PCA",
                             "MI_CART", "MI_RF", "mean")
    )
  } )
  bias_plots_c2 <- lapply(1:length(bias_res_list), function(p){
    plot_gg(bias_res_list[[p]][[2]], 
            type = "bias_worst",
            plot_name = c("High Dimensional Condition"),
            meth_compare = c("DURR_la", "IURR_la", "blasso", "bridge",
                             "MI_PCA",
                             "MI_CART", "MI_RF", "mean")
    )
  } )
  
  # Display plots
  p1 <- arrangeGrob(grobs = list(bias_plots_c1[[1]], bias_plots_c2[[1]]),
                    top = "Model 1: \u03B2 Religiosity", ncol = 2,
  )
  p2 <- arrangeGrob(grobs = list(bias_plots_c1[[2]], bias_plots_c2[[2]]),
                    top = "Model 2: \u03B2 Nativist Attitudes", ncol = 2
  )
  pf <- grid.arrange(p1, p2)
  ggsave(file = "../output/graphs/exp4_imp_bias_worst.pdf", pf)
  
# > CI   ####
  
  # Extract what you want to plot
  ci_res_list <- list(lapply(1:length(res$m1),
                               function(x) res$m1[[x]]$ci_cov["rel", ]),
                        lapply(1:length(res$m2),
                               function(x) res$m2[[x]]$ci_cov["NatAt", ]))
  
  # Dfine plots
  ci_plots_c1 <- lapply(1:length(ci_res_list), function(p){
    plot_gg(dt        = ci_res_list[[p]][[1]], 
            type      = "ci",
            plot_name = c("Low Dimensional Condition"),
            meth_compare = c("DURR_la", "IURR_la", "blasso", "bridge",
                             "MI_PCA",
                             "MI_CART", "MI_RF", "mean")
    ) } )
  ci_plots_c2 <- lapply(1:length(ci_res_list), function(p){
    plot_gg(dt        = ci_res_list[[p]][[2]], 
            type      = "ci",
            plot_name = c("High Dimensional Condition"),
            meth_compare = c("DURR_la", "IURR_la", "blasso", "bridge",
                             "MI_PCA",
                             "MI_CART", "MI_RF", "mean")
    ) } )
  
  # Combine and arrange plots
  p1 <- arrangeGrob(grobs = list(ci_plots_c1[[1]], ci_plots_c2[[1]]),
                     top = "Model 1: \u03B2 Religiosity", ncol = 2
  )
  
  p2 <- arrangeGrob(grobs = list(ci_plots_c1[[2]], ci_plots_c2[[2]]),
                    top = "Confidence Interval Coverage (CIC) Nativist Attitudes' \u03B2",
                    ncol = 2
  )
  pf <- grid.arrange(p1, p2)
  pf <- grid.arrange(p2) # single plot
  ggsave(file = "../output/graphs/exp4_imp_ci.pdf", pf, device = cairo_pdf)
  
  ## Worst CI coverage
  # Extract what you want to plot
  ci_res_list <- list(lapply(1:length(res$m1),
                             function(x) res$m1[[x]]$ci_cov),
                      lapply(1:length(res$m2),
                             function(x) res$m2[[x]]$ci_cov))
  
  # Dfine plots
  ci_plots_c1 <- lapply(1:length(bias_res_list), function(p){
    plot_gg(dt        = ci_res_list[[p]][[1]], 
            type      = "ci_worst",
            plot_name = c("Low Dimensional Condition"),
            meth_compare = c("DURR_la", "IURR_la", "blasso", "bridge",
                             "MI_PCA",
                             "MI_CART", "MI_RF", "mean")
    ) } )
  ci_plots_c2 <- lapply(1:length(bias_res_list), function(p){
    plot_gg(dt        = ci_res_list[[p]][[2]], 
            type      = "ci_worst",
            plot_name = c("High Dimensional Condition"),
            meth_compare = c("DURR_la", "IURR_la", "blasso", "bridge",
                             "MI_PCA",
                             "MI_CART", "MI_RF", "mean")
    ) } )
  
  # Combine and arrange plots
  p1 <- arrangeGrob(grobs = list(ci_plots_c1[[1]], ci_plots_c2[[1]]),
                    top = "Model 1: \u03B2 Religiosity", ncol = 2
  )
  
  p2 <- arrangeGrob(grobs = list(ci_plots_c1[[2]], ci_plots_c2[[2]]),
                    top = "Model 2: \u03B2 Nativist Attitudes", ncol = 2
  )
  pf <- grid.arrange(p1, p2)
  ggsave(file = "../output/graphs/exp4_imp_ci_worst.pdf", pf)
  
  ## > CI Width ####
  ci_res_list <- list(sapply(1:length(res$m1),
                             function(x) res$m1[[x]]$CIW[9, ]),
                      sapply(1:length(res$m2),
                             function(x) res$m2[[x]]$CIW[3, ]))
  
# Multivariate distance ---------------------------------------------------
# > Bias ####
  # Per model
  # Extract what you want to plot
  bias_res_list <- list(lapply(1:length(res$m1),
                               function(x) res$m1[[x]]$ed_est),
                        lapply(1:length(res$m2),
                               function(x) res$m2[[x]]$ed_est))
  
  # Define plots
  bias_plots_c1 <- lapply(1:length(bias_res_list), function(p){
    plot_gg(dt        = bias_res_list[[p]][[1]], 
            type      = "ed",
            plot_name = c("Low Dimensional Condition"),
            meth_compare = c("DURR_la", "IURR_la", "blasso", "bridge",
                             "MI_PCA",
                             "MI_CART", "MI_RF", "mean", "CC")
    ) } )
  bias_plots_c2 <- lapply(1:length(bias_res_list), function(p){
    plot_gg(dt        = bias_res_list[[p]][[2]], 
            type      = "ed",
            plot_name = c("High Dimensional Condition"),
            meth_compare = c("DURR_la", "IURR_la", "blasso", "bridge",
                             "MI_PCA",
                             "MI_CART", "MI_RF", "mean", "CC")
    ) } )
  
  # Combine and arrange plots
  p1 <- arrangeGrob(grobs = list(bias_plots_c1[[1]], bias_plots_c2[[1]]),
                     top = "Model 1", ncol = 2
  )
  p2 <- arrangeGrob(grobs = list(bias_plots_c1[[2]], bias_plots_c2[[2]]),
                     top = "Euclidean Distance between vector of true and after imputation estimates", 
                    ncol = 2
  )
  pf <- grid.arrange(p1, p2)
  pf <- grid.arrange(p2) # single plot
  ggsave(file = "../output/graphs/exp4_ed_bias.pdf", pf)
  
  ## OVERALL ##
  # Extract what you want to plot
  bias_res_list <- list(res$ed_all[[1]]$ed_est,
                        res$ed_all[[2]]$ed_est)
  
  # Define plots
  bias_plots_c1 <- plot_gg(dt        = bias_res_list[[1]], 
                           type      = "ed",
                           plot_name = c("(a)"),
                           meth_compare = c("DURR_la", "IURR_la", "blasso", "bridge",
                                            "MI_PCA",
                                            "MI_CART", "MI_RF", "mean")
  )
  bias_plots_c2 <- plot_gg(dt        = bias_res_list[[2]], 
                           type      = "ed",
                           plot_name = c("(b)"),
                           meth_compare = c("DURR_la", "IURR_la", "blasso", "bridge",
                                            "MI_PCA",
                                            "MI_CART", "MI_RF", "mean")
  )
  pf <- grid.arrange(bias_plots_c1, bias_plots_c2, ncol = 2)
  # ggsave(file = "../output/graphs/exp4_ed_bias.pdf", pf)
  
# > CI ####
  # Extract what you want to plot
  ci_res_list <- list(lapply(1:length(res$m1),
                             function(x) res$m1[[x]]$ed_ci),
                      lapply(1:length(res$m2),
                             function(x) res$m2[[x]]$ed_ci))
  
  # Dfine plots
  ci_plots_c1 <- lapply(1:length(ci_res_list), function(p){
    plot_gg(dt        = ci_res_list[[p]][[1]], 
            type      = "ed",
            plot_name = c("Low Dimensional Condition"),
            meth_compare = c("DURR_la", "IURR_la", "blasso", "bridge",
                             "MI_PCA",
                             "MI_CART", "MI_RF", "mean")
    ) } )
  ci_plots_c2 <- lapply(1:length(ci_res_list), function(p){
    plot_gg(dt        = ci_res_list[[p]][[2]], 
            type      = "ed",
            plot_name = c("High Dimensional Condition"),
            meth_compare = c("DURR_la", "IURR_la", "blasso", "bridge",
                             "MI_PCA",
                             "MI_CART", "MI_RF", "mean")
    ) } )
  
  # Combine and arrange plots
  p1 <- arrangeGrob(grobs = list(ci_plots_c1[[1]], ci_plots_c2[[1]]),
                     top = "Model 1", ncol = 2
  )
  p2 <- arrangeGrob(grobs = list(ci_plots_c1[[2]], ci_plots_c2[[2]]),
                     top = "Model 2", ncol = 2
  )
  pf <- grid.arrange(p1, p2)
  pf <- grid.arrange(p2) # single plot
  ggsave(file = "../output/graphs/exp4_ed_ci.pdf", pf)
  
  ## OVERALL ##
  # Extract what you want to plot
  ci_res_list <- list(res$ed_all[[1]]$ed_ci,
                      res$ed_all[[2]]$ed_ci)
  
  # Define plots
  ci_plots_c1 <- plot_gg(dt        = ci_res_list[[1]], 
                         type      = "ed",
                         plot_name = c("(a)"),
                         meth_compare = c("DURR_la", "IURR_la", "blasso", "bridge",
                                          "MI_PCA",
                                          "MI_CART", "MI_RF", "mean")
  )
  ci_plots_c2 <- plot_gg(dt        = ci_res_list[[2]], 
                         type      = "ed",
                         plot_name = c("(b)"),
                         meth_compare = c("DURR_la", "IURR_la", "blasso", "bridge",
                                          "MI_PCA",
                                          "MI_CART", "MI_RF", "mean")
  )
  pf <- grid.arrange(ci_plots_c1, ci_plots_c2, ncol = 2)
  # ggsave(file = "../output/graphs/exp4_ed_bias.pdf", pf)
  
  
  ## > Confidence Interval Width ####
  ci_res_list <- list(lapply(1:length(res$m1),
                             function(x) res$m1[[x]]$CIW),
                      lapply(1:length(res$m2),
                             function(x) res$m2[[x]]$CIW))
  lapply(1:length(ci_res_list), function(i){
    t(sapply(1:length(ci_res_list[[i]]), function(j){
      # round(colMeans(ci_res_list[[i]][[j]]), 3)
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
  # Read R object
  out <- readRDS(paste0("../output/", filename, ".rds"))
  
  # Produce data for plot
  out_time <- sapply(1:length(names(out[[1]])), res_sem_time, out = out)
  colnames(out_time) <- names(out[[1]])
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
    theme(plot.title = element_text(hjust = 0.5),
          axis.text = element_text(size = 12)) +
    
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
               # linetype = "longdash", 
               color = "gray") +
    coord_cartesian(xlim = plot_xlim)
  
  return(p)
  
  }
  
  ## Use fucntion on time data
  p <- lapply(1:ncol(out_time), function(j){
    plot_gg_time( round(out_time[, j],1),
                  plot_name = c("Low Dimensional Condition", 
                                "High Dimensional Condition")[j],
                  meth_compare = c("DURR_la", "IURR_la", "blasso", "bridge",
                                   "MI_PCA",
                                   "MI_CART", "MI_RF")
                  )
  })
  
  p <- arrangeGrob(grobs = list(p[[1]], p[[2]]),
                   top = "Average Imputation time (minutes) for each MI method",
                   ncol = 2)
  pf <- grid.arrange(p)

  ggsave(file = "../output/graphs/exp4_time.pdf", pf)

