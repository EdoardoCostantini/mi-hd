res <- output
res$conds
# Full information names
# cond_parms <- paste0("l = ",  res$conds$lv, ", ",
#                      "pm = ",  res$conds$pm)
cond_parms <- paste0("l = ",  res$conds$lv, ", ",
                     "pm = ",  res$conds$pm, "\n",
                     "\u03bb = ", c(rep("(0.9, 0.97)", 4)))

cond_labels <- c("low-dim-low-pm-high-\u03bb",
                 "high-dim-low-pm-high-\u03bb",
                 "low-dim-high-pm-high-\u03bb",
                 "high-dim-high-pm-high-\u03bb")
                 # "low-dim-low-pm-low-\u03bb",
                 # "high-dim-low-pm-low-\u03bb",
                 # "low-dim-high-pm-low-\u03bb",
                 # "high-dim-high-pm-low-\u03bb")
cond_names <- paste(cond_labels, cond_parms, sep = " \n ")
# Numbered Names
# cond_names <- paste0("Condition ", 1:8)
data.frame(res$conds, cond_names)

# Select conditions to print
conds_select <- 5:8
conds_select <- 1:4

# Bias (Facet grid) -------------------------------------------------------

meth_compare = c("DURR_all","DURR_si","IURR_all","IURR_si","blasso",
                 #"bridge",
                 "MI_PCA","MI_CART" ,"MI_RF","MI_OP",
                 "missFor","CC")
meth_compare = c("blasso", "MI_OP", "mean", "CC", "GS")

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