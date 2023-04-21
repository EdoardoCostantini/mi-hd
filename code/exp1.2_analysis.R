# Project:   imputeHD-comp
# Objective: Analysis of results for experiment 1.2
# Author:    Edoardo Costantini
# Created:   2023-04-12
# Modified:  2023-04-21
# Notes: 

rm(list = ls())
source("./init_general.R")

# Read results from the combined results
res <- readRDS("../output/exp1_2_simOut_20230415_0830_res.rds")
res <- readRDS("../output/exp1_2_simOut_20230419_1102_res.rds")
res <- readRDS("../output/exp1_2_simOut_20230419_1510_res.rds")
res <- readRDS("../output/exp1_2_simOut_20230419_1403_res.rds")

# Change names of methods if required
levels(res$methods) <- str_replace(levels(res$methods), "blasso", "BLasso")
levels(res$methods) <- str_replace(levels(res$methods), "bridge", "BRidge")
levels(res$methods) <- str_replace(levels(res$methods), "MI-qp", "MI-QP")
levels(res$methods) <- str_replace(levels(res$methods), "MI-am", "MI-AM")
levels(res$methods) <- str_replace(levels(res$methods), "MI-OP", "MI-OR")
levels(res$methods) <- str_replace(levels(res$methods), "stepFor", "MI-SF")

# Which plot to plot
p_grep <- c("50", "500")[2]

# Plot Sizes Decisions

gp_width <- 15
gp_height <- 20
sp_width <- 3.75
sp_height <- 6

# Plot font

plot_text_family <- "sans"
plot_text_face <- "plain"
plot_text_size <- 9

# 

output_sem <- res

# Condition Names / Labels
label_cond <- unique(output_sem$cond)
label_parm <- unique(output_sem$parm)

# Threshold lines
vline_bias <- 10

# SE for threshold
ci_lvl <- .95
dt_reps <- 1e3
SEp <- sqrt(ci_lvl * (1 - ci_lvl) / dt_reps)
low_thr <- (.95 - SEp * 2) * 100
hig_thr <- (.95 + SEp * 2) * 100
vline_burton <- c(low_thr, hig_thr)
vline_vanBuu <- 90
xci_breaks <- sort(c(0, 10, 20, 50))

# Fix methods order
# output_sem$methods <- factor(output_sem$methods,
#     levels = levels(output_sem$methods)[c(1:7, 13, 11:12, 8, 9, 10)]
# )

# Bias
x <- 1 # bias
methods_sel <- levels(output_sem$methods)[1:11]

pf <- output_sem %>%
    filter(
        analysis == unique(analysis)[x],
        grepl(p_grep, cond),
        variable %in% c("Min", "Mean", "Max"),
        methods %in% methods_sel
    ) %>%
    mutate(methods = fct_relabel(
        methods,
        str_replace,
        "-la", ""
    )) %>%
    # Main Plot
    ggplot(data = ., aes(
        y = methods,
        x = value,
        shape = variable
    )) +
    geom_point(size = 1.75) +
    geom_line(aes(group = methods),
        size = .25
    ) +
    # Grid
    facet_grid(
        rows = vars(factor(parm,
            levels = unique(parm)
        )),
        cols = vars(cond)
    ) +
    geom_vline(
        data = data.frame(
            xint = 10,
            analysis = "Percentage Relative Bias"
        ),
        linetype = "solid",
        size = .15,
        aes(xintercept = xint)
    ) +

    # Format
    scale_x_continuous(
        labels = xci_breaks,
        breaks = xci_breaks
    ) +
    scale_y_discrete(limits = rev) +
    scale_shape_manual(values = c("I", "I", "I")) +
    coord_cartesian(xlim = c(0, 50)) +
    labs( # title = label_parm[x],
        x = NULL,
        y = NULL,
        linetype = NULL,
        shape = NULL
    ) +
    theme(
        panel.background = element_rect(
            fill = NA,
            color = "gray"
        ),
        panel.grid.major = element_line(
            color = "gray",
            size = 0.15,
            linetype = 1
        ),
        legend.key = element_rect(
            colour = "gray",
            fill = NA,
            size = .15
        ),
        text = element_text(
            family = plot_text_family,
            face = plot_text_face,
            size = plot_text_size
        ),
        axis.ticks = element_blank(),
        legend.position = "none"
    )

pf

plot_name <- paste0(
    "../output/graphs/exp1.2_", "bias", "_",
    p_grep, "revision.pdf"
)
ggsave(
    file = plot_name,
    width = gp_width, height = gp_height * .65,
    units = "cm",
    pf
)

# Confidence interval coverage ------------------------------------------------

x <- 2 # Confidence intervals
methods_sel <- levels(output_sem$methods) # [-8]

# SE for threshold
ci_lvl <- .95
dt_reps <- 500
SEp <- sqrt(ci_lvl * (1 - ci_lvl) / dt_reps)
low_thr <- (.95 - SEp * 2) * 100
hig_thr <- (.95 + SEp * 2) * 100
vline_burton <- c(low_thr, hig_thr)
vline_vanBuu <- c(90, 99)
xci_breaks <- sort(c(80, vline_vanBuu, 95, round(vline_burton, 1), 100))

pf <- output_sem %>%
  filter(
    analysis == unique(analysis)[x],
    grepl(p_grep, cond),
    variable %in% c("Min", "Mean", "Max"),
    methods %in% methods_sel
  ) %>%
  # Drop pm = ** as we are plotting only one value for this condition
  mutate(cond = fct_relabel(
    cond,
    str_replace,
    " pm = [0-9]\\.[0-9]", ""
  )) %>%
  mutate(methods = fct_relabel(
    methods,
    str_replace,
    "-la", ""
  )) %>%
  # Main Plot
  ggplot(data = ., aes(
    y = methods,
    x = value,
    shape = variable
  )) +
  geom_point(size = 1.75, show.legend = FALSE) +
  geom_line(aes(group = methods),
    size = .25
  ) +
  # Grid
  facet_grid(
    rows = vars(factor(parm,
      levels = unique(parm)
    )),
    cols = vars(cond)
  ) +
  geom_vline(
    data = data.frame(
      xint = 95,
      analysis = "CI coverage"
    ),
    linetype = "solid",
    size = .15,
    aes(
      xintercept = xint,
      lty = paste0("nominal level")
    )
  ) + 

  # Format
  scale_y_discrete(limits = rev) +
  scale_x_continuous(
    labels = as.character(round(xci_breaks / 100, 2)),
    breaks = xci_breaks
  ) +
  coord_cartesian(xlim = c(min(xci_breaks), max(xci_breaks))) +
  scale_shape_manual(values = c("I", "I", "I")) +
  # scale_linetype_manual(values = c("longdash", "dashed", "solid")) +
  # guides(linetype = guide_legend(override.aes = list(size = c(.20, .20, .15)))) +
  labs( # title = label_parm[x],
    x = NULL,
    y = NULL,
    linetype = NULL,
    shape = NULL
  ) +
  theme(
    panel.background = element_rect(
      fill = NA,
      color = "gray"
    ),
    panel.grid.major = element_line(
      color = "gray",
      size = 0.175,
      linetype = 1
    ),
    axis.ticks = element_blank(),
    legend.key = element_rect(
      colour = "gray",
      fill = NA,
      size = .15
    ),
    legend.position = "bottom",
    text = element_text(
      family = plot_text_family,
      face = plot_text_face,
      size = plot_text_size
    ),
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    )
  )

pf

plot_name <- paste0(
  "../output/graphs/exp1.2_", "CI", "_",
  p_grep, "revision.pdf"
)
ggsave(
  file = plot_name,
  width = gp_width, height = gp_height * .65,
  units = "cm",
  pf
)

# Confidence Interval Width ---------------------------------------------------

# SE for threshold
dt_reps <- 500

# Bias
x <- 3 # CIW
methods_sel <- levels(output_sem$methods)[-c(12)]

pf <- output_sem %>%
    filter(
        analysis == unique(analysis)[x],
        grepl(p_grep, cond),
        variable %in% c("Mean"),
        methods %in% methods_sel
    ) %>%
    # Drop pm = ** as we are plotting only one value for this condition
    mutate(cond = fct_relabel(
        cond,
        str_replace,
        " pm = [0-9]\\.[0-9]", ""
    )) %>%
    mutate(methods = fct_relabel(
        methods,
        str_replace,
        "-la", ""
    )) %>%
    # Main Plot
    ggplot(data = ., aes(
        y = methods,
        x = value
    )) +
    geom_point(size = 1.75) +
    # Grid
    facet_grid(
        rows = vars(factor(parm,
            levels = unique(parm)
        )),
        cols = vars(cond)
    ) +
    coord_cartesian(xlim = c(0, 10)) +

    # Format
    scale_y_discrete(limits = rev) +
    labs( # title = label_parm[x],
        x = NULL,
        y = NULL,
        linetype = NULL
    ) +
    theme(
        panel.background = element_rect(
            fill = NA,
            color = "gray"
        ),
        panel.grid.major = element_line(
            color = "gray",
            size = 0.15,
            linetype = 1
        ),
        legend.key = element_rect(
            colour = "gray",
            fill = NA,
            size = .15
        ),
        text = element_text(
            family = plot_text_family,
            face = plot_text_face,
            size = plot_text_size
        ),
        axis.ticks = element_blank(),
        legend.position = "none"
    )

pf

plot_name <- paste0(
    "../output/graphs/exp1.2_", "ciw", "_",
    p_grep, "revision.pdf"
)

ggsave(
    file = plot_name,
    width = gp_width, height = gp_height * .65,
    units = "cm",
    pf
)

# Time plot ---------------------------------------------------------------

dt_outtime <- readRDS("../output/exp1_2_simOut_20230408_1748_time.rds")

round(dt_outtime, 1)
round(dt_outtime/60, 0)

