### Title:    Analysis of results from experiment 2
### Project:  Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-08-27
### Notes:    reads output form results.R script and shows the numbers that
###           are used to draw the conclusions.

  rm(list = ls())
  source("./init_general.R")

# Read results from a run of simulation study
  # Result selection
  filename <- "exp4_simOut_20201016_2341" # old
  filename <- "exp4_simOut_20201019_1344" # current one
  filename <- "exp4_simOut_20201027_1610" # combined results from more runs
  filename <- "exp4_simOut_20201204_2121" # updated model 1, 500 data
  
  # Read R object
  res <- readRDS(paste0("../output/", filename, "_res.rds"))
  out <- res
  
  # Plot Size
  sp_width <- 5
  sp_height <- 4
  
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
  
# Variable of interest ----------------------------------------------------

  # > Bias (Facet) ####
  # Which Methods do you want to plot
  meths <- c("DURR_la", "IURR_la", 
             "blasso", "bridge",
             "MI_PCA",
             "MI_CART", "MI_RF", 
             "missFor", 
             "mean", "CC")
  pf <- plot_exp4(dt = list(lapply(1:length(res$m1),
                                   function(x) res$m1[[x]]$bias_per["rel", ]),
                            lapply(1:length(res$m2),
                                   function(x) res$m2[[x]]$bias_per["NatAt", ])),
            type = "bias",
            dt_reps = 500,
            ci_lvl = .95,
            plot_cond = NULL,
            plot_name = NULL,
            bar_col = "#595959",
            meth_compare = meths,
            meth_sort = FALSE)
  pf
  
  ggsave(pf,
         file = "../output/graphs/exp4_imp_bias.pdf",
         width = sp_width*2, height = sp_height*2,
         units = "cm",
         device = cairo_pdf)
  
# > Bias ####

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
         width = 8.25/3*2, height = 11.75/4*2,
         device = cairo_pdf)
  
  
  # > CI (Facet) ####
  # Which Methods do you want to plot
  meths <- c("DURR_la", "IURR_la", 
             "blasso", "bridge",
             "MI_PCA",
             "MI_CART", "MI_RF", 
             "missFor", 
             "mean", "CC", "GS")
  
  pf <- plot_exp4(dt = list(lapply(1:length(res$m1),
                             function(x) res$m1[[x]]$ci_cov["rel", ]),
                      lapply(1:length(res$m2),
                             function(x) res$m2[[x]]$ci_cov["NatAt", ])),
            type = "ci",
            dt_reps = 500,
            ci_lvl = .95,
            plot_cond = NULL,
            plot_name = NULL,
            bar_col = "#595959",
            meth_compare = meths,
            meth_sort = FALSE)
  ggsave(pf, 
         file = "../output/graphs/exp4_imp_ci.pdf", 
         width = sp_width*2, height = sp_height*2,
         units = "cm",
         device = cairo_pdf)
  
# > CI   ####
  # Which Methods do you want to plot
  meths <- c("DURR_la", "IURR_la", 
             "blasso", "bridge",
             "MI_PCA",
             "MI_CART", "MI_RF", 
             "missFor", 
             "mean", "CC", "GS")
  
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
            plot_name = "",
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
         width = 8.25/3*2, height = 11.75/4*2,
         device = cairo_pdf)
  

  # > CIW (Facet) ####
  # Which Methods do you want to plot
  meths <- c("DURR_la", "IURR_la", 
             "blasso", "bridge",
             "MI_PCA",
             "MI_CART", "MI_RF", 
             "missFor", 
             "mean", "CC", "GS")
  
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
  
  # > CIW   ####
  # Define plots
  ci_plots_m1 <- lapply(1:2, function(p){
    plot_gg(
            # dt        = lapply(1:length(res$m1),
            #                    function(x) res$m1[[x]]$CIW["rel", ])[[p]],
            dt        = data.frame(t(colMeans(lapply(1:length(res$m1),
                                                     function(x) res$m1[[x]]$CIW)[[p]]))),
            type      = "ciw",
            plot_name = cond_names[[p]],
            meth_compare = meths) } )
  
  ci_plots_m2 <- lapply(1:2, function(p){
    plot_gg(
            # dt        = lapply(1:length(res$m2),
            #                    function(x) res$m2[[x]]$CIW["NatAt", ])[[p]],
            dt        = data.frame(t(colMeans(lapply(1:length(res$m2),
                                                     function(x) res$m2[[x]]$CIW)[[p]]))),
            type      = "ciw",
            plot_name = "",
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
         width = 8.25/3*2, height = 11.75/4*2,
         device = cairo_pdf)
  
# Multivariate distance ---------------------------------------------------
  
  # > ED BIAS (Facet) ####
  # Which Methods do you want to plot
  meths <- c("DURR_la", "IURR_la", 
             "blasso", "bridge",
             "MI_PCA",
             "MI_CART", "MI_RF", 
             "missFor", 
             "mean", "CC")
  
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
  
# > Bias ####
  meths <- c("DURR_la", "IURR_la", 
             "blasso", "bridge",
             "MI_PCA",
             "MI_CART", "MI_RF", 
             "missFor", 
             "mean", "CC")
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
            plot_name = "",
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
  # grid.arrange(p2) # single plot
  
  # Save Plot
  ggsave(pf,
         file = "../output/graphs/exp4_ed_bias.pdf", 
         width = 8.25/3*2, height = 11.75/4*2)

# > ED CI (Facet) ####
  # Which Methods do you want to plot
  meths <- c("DURR_la", "IURR_la", 
             "blasso", "bridge",
             "MI_PCA",
             "MI_CART", "MI_RF", 
             "missFor", 
             "mean", "CC", "GS")
  
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
  ggsave(file = "../output/graphs/exp4_ed_ci.pdf", pf,
         width = sp_width*2, height = sp_height*2,
         units = "cm")
  
# > CI ####
  meths <- c("DURR_la", "IURR_la", 
             "blasso", "bridge",
             "MI_PCA",
             "MI_CART", "MI_RF", 
             "missFor", 
             "mean", "CC", "GS")
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
            plot_name = "",
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
         width = 8.25/3*2, height = 11.75/4*2)

# Time plot ---------------------------------------------------------------

  rm(list = ls())
  source("./init_general.R")
  
  # Result selection
  filename <- "exp4_simOut_20201019_1344" # current one
  filename <- "exp4_simOut_20201027_1610" # current one
  filename <- "exp4_simOut_20201204_2121" # updated model 1, 500 data
  
  # Read R object
  res <- readRDS(paste0("../output/", filename, ".rds"))

  # Plot Size
  sp_width <- 5
  sp_height <- 4
  
  # Produce data for plot
  out_time <- sapply(1:length(names(res[[1]])), res_sem_time, out = res)
  colnames(out_time) <- names(res[[1]])
  dt = t(out_time)
  
  # Get Plot
  pf <- plot_time(dt = t(out_time), 
            plot_cond = NULL,
            plot_name = NULL,
            meth_compare = rev(c("DURR_la", "IURR_la", "blasso", "bridge",
                                 "MI_PCA",
                                 "MI_CART", "MI_RF")),
            meth_sort = FALSE)
  pf
  # Save Plot  
  ggsave(pf,
         file = "../output/graphs/exp4_time.pdf", 
         width = sp_width*2, height = sp_height*1,
         units = "cm")
