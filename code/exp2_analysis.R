### Title:    Analysis of results from experiment 2
### Project:  Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-07-09
### Notes:    reads output form results.R script and shows the numbers that
###           are used to draw the conclusions.

    library(xtable)
    
    rm(list = ls())
    source("./init_general.R")
  
  # Read results from a run of simulation study
    # exp2_res <- readRDS("../output/exp2_simOut_20200812_1449_res.rds") # way out
    filename <- "exp2_simOut_20200819_1743"
    out <- readRDS(paste0("../output/", filename, ".rds"))
    exp2_res <- readRDS(paste0("../output/", filename, "_res.rds"))

# Bias --------------------------------------------------------------------
  
# Recap of set up  

  out$conds
  
  # c  lv  pm   fl ridge
  # 1  10 0.1 high 1e-01
  # 2 100 0.1 high 1e-07
  # 3  10 0.3 high 1e-01
  # 4 100 0.3 high 1e-07
  # ------------------ #
  # 5  10 0.1  low 1e-01
  # 6 100 0.1  low 1e-07
  # 7  10 0.3  low 1e-01
  # 8 100 0.3  low 1e-07
  
  list(
    it_number = out$parms$n_it,
    mis_var   = out$parms$z_m_id,
    miss_type = out$parms$missType,
    lv_rm_x   = out$parms$rm_x,
    lv_number = "condition specific: 10 or 100",
    factor_loadings = "a) runif btw .9 and .97; b) runif btw .5 and .6",
    substantive = "SEM, CFA raw data; SEM, LM scored data (mean of items)"
  )
  
  # Condition indexes
  cindex_lh <- c(1:4)
  cindex_hd <- c(2, 4)
  cindex_hp <- c(3, 4)
  
#> SEM scored ####
  # PCA issue with variances remains but good performances
  # HIGH FACTOR LOADINGS
  lapply(exp2_res$semS,
         function(x) x$bias_sd)[cindex_lh] # means in terms of sd bias
  lapply(exp2_res$semS,
         function(x) x$bias_per[-c(1:2),])[cindex_lh] # Var Covar percent bias
  # LOW FACTOR LOADINGS
  lapply(exp2_res$semS,
         function(x) x$bias_sd)[-cindex_lh]
  lapply(exp2_res$semS,
         function(x) x$bias_per[-c(1:2),])[-cindex_lh]
  
  # CONFIDENCE INTERVALS
  # Look at sizes
  lapply(exp2_res$semS,
         function(x) x$ci_cov)[c(4,8)]
  # Look at ED
  table.ed.obj <- round(t(sapply(exp2_res$semS,
                           res_ed_ci)), 0)
  
  # Table
  table.ed <- xtable(table.ed.obj, 
                     align = c("l",# 2nd slot is for actual rownames
                               rep("c", ncol(table.ed.obj)) ),
                     digits = 0,
                     caption = "Euclidean distance of the vector of 95\\% CIs"
  )
  print(table.ed, 
        hline.after = c(0, 4,
                        nrow(table.ed.obj)), # repeating horizonta
        # add.to.row = addtorow, 
        include.rownames = TRUE)
  
  # Fix Column names and the like
  addtorow         <- list()
  addtorow$pos     <- list(0, 0, 4, 4) # Add Condition 2 break
  addtorow$command <- c(
    # Add line break
    paste0('\\hline \\\\\n'),
    # factor loadings separation
    paste0(paste0('& \\multicolumn{',
                  ncol(table.ed)-1,
                  '}{c}{high $\\lambda$}', 
                  collapse=''), 
           '\\\\\n'),
    # Add line break
    paste0('\\hline \\\\\n'),
    # Condition Break
    paste0(paste0('& \\multicolumn{',
                  ncol(table.ed)-1,
                  '}{c}{low $\\lambda$}', 
                  collapse=''), 
           '\\\\\n')
  )
  
  print(table.ed, 
        hline.after = c(0, 4,
                        nrow(table.ed.obj)), # repeating horizonta
        add.to.row = addtorow,
        include.rownames = TRUE,
        scalebox = 0.5)
  
#> SEM raw data ####
  # t(sapply(exp2_res$semR,
  #          function(x) x$validReps))
  
  # EUCLIDEAN distance measure for means
  # LOW FACTOR LOADINGS
  t(sapply(exp2_res$semR,
           res_ed_est, index = 1:10))[5:8, ]
  
  # VARIANCES
  # Good PCA good for all but highest condition with low factor loadings
  # Last condition is the interesting one: PCA and IURR perform well in all
  # other conditions but when factor loadings are smaller then PCA starts to
  # show its biased variances problem, and IURR its biased covariances problem
  lapply(exp2_res$semR,
         function(x) x$bias_per[-c(1:10),])[c(4,8)]
  
  # PCA variance struggle
  t(sapply(exp2_res$semR,
           res_ed_est, index = 11:20))
  
  # PCA covariances domination
  t(sapply(exp2_res$semR,
           res_ed_est, index = -c(1:20)))
  
  # Confidence Intervals one by one too crowded
  # lapply(exp2_res$semR,
  #        function(x) x$ci_cov)[c(4, 8)]
  
  # Look at Euclidean distance measure
  # CI are great for PCA!
  t(sapply(exp2_res$semR,
           res_ed_ci))[cindex_lh, ]
  t(sapply(exp2_res$semR,
           res_ed_ci))[-cindex_lh, ]
  
#> CFA raw data ####
  # FACTOR LOADINGS: PCA does not even budge with high dim and low factor 
  # loadings
  lapply(exp2_res$CFA,
         function(x) x$bias_per[1:10, ])[cindex_lh]
  lapply(exp2_res$CFA,
         function(x) x$bias_per[1:10, ])[-cindex_lh]
  # # CONFIDENCE INTERVALS ARE GOOD AS WELL
  # lapply(exp2_res$CFA,
  #        function(x) x$ci_cov[1:10, ])[cindex_lh]
  # lapply(exp2_res$CFA,
  #        function(x) x$ci_cov[1:10, ])[-cindex_lh]
  
  # Latex table
  lapply(exp2_res$CFA,
         function(x) x$bias_per[1:10, ])[1]
  
  # Same pm, different lv
  cond.seq <- 3:4
  parm.seq <- 1:5 # for first latent variable
  meth.seq <- 1:8
  conds <- out$conds
  # Structure results for MI methods for selected conditions, parameters, and methods
  mi_bici <- lapply(cond.seq, function(x){
    cond.obj   <- conds[x, ]
    target.obj <- exp2_res$CFA[[x]]
    temp <- rbind(Bias = t(target.obj$bias_per[parm.seq, 
                                               meth.seq+1]),
                  cov = t(round(target.obj$ci_cov[parm.seq, 
                                                  meth.seq], 0)))
    meth.names <- unique(rownames(temp))
    temp <- structure(temp, dim = c(length(meth.seq), length(parm.seq)*2))
    # temp.df <- data.frame(temp)
    # rownames(temp.df) <- meth.names
    mi_bici.out <- cbind(meths = meth.names, temp)
    return(mi_bici.out)
  })
  mi_bici <- do.call(rbind, mi_bici)
  
  # Structure results for SI methods for selected conditions, parameters, and methods
  meth.seq <- 9:11
  si_bici <- lapply(cond.seq[1], function(x){
    cond.obj   <- conds[x, ]
    target.obj <- exp2_res$CFA[[x]]
    temp <- rbind(Bias = t(target.obj$bias_per[parm.seq, 
                                               meth.seq+1]),
                  cov = t(round(target.obj$ci_cov[parm.seq, 
                                                  meth.seq], 0)))
    meth.names <- unique(rownames(temp))
    temp <- structure(temp, dim = c(length(meth.seq), length(parm.seq)*2))
    # temp.df <- data.frame(temp)
    # rownames(temp.df) <- meth.names
    si_bici.out <- cbind(meths = meth.names, temp)
    return(si_bici.out)
  })
    
  si_bici <- do.call(rbind, si_bici)

  # Create base table
  table.1 <- xtable(rbind(si_bici, mi_bici), 
                    align = c("", # 1st slot is for R object rownames (not displayed)
                              "l",# 2nd slot is for actual rownames
                              rep("c", ncol(mi_bici)-1) )
                    )
  
  # Fix Column names and the like
  addtorow         <- list()
  addtorow$pos     <- list(0, # Add latent factors header
                           0, # Add type of information reported
                           3, # line break
                           3, # Add Condition 1 break
                           11, # line break
                           11) # Add Condition 2 break
  addtorow$command <- c(
    # Factor loading name
    paste0(paste0('& \\multicolumn{2}{c}{', 
           paste0("$\\lambda_", parm.seq, "$"), 
           '}', 
           collapse=''), '\\\\\n'),
    # Bias / Cov column names
    paste0(
      paste0(" & ", rep(c("bias", "cov"), (ncol(table.1)-1)/2), collapse = ""), 
           '\\\\\n'),
    # Add line break
    paste0('\\hline \\\\\n'),
    # Condition Break
    paste0(paste0('& \\multicolumn{',
                  ncol(table.1)-1,
                  '}{c}{', 
                  paste0("pm = ", conds[cond.seq[1], ]$pm,
                         ", lv = ", conds[cond.seq[1], ]$lv),
                  '}', 
                  collapse=''), 
           '\\\\\n'),
    # Add line break
    paste0('\\hline \\\\\n'),
    # Condition Break 
    paste0(paste0('& \\multicolumn{',
                  ncol(table.1)-1,
                  '}{c}{', 
                  paste0("pm = ", conds[cond.seq[2], ]$pm,
                         ", lv = ", conds[cond.seq[2], ]$lv),
                  '}', 
                  collapse=''), 
           '\\\\\n')
  )
  
  print(table.1, 
        hline.after = c(0, 
                        3,
                        11), # repeating horizonta
        add.to.row = addtorow, 
        include.colnames = FALSE,
        include.rownames = FALSE,
        scalebox = 0.5)
  
  print.xtable()
  # 
  # temp <- rbind(Bias = t(exp2_res$CFA$cond1$bias_per[1:3, 2:3]),
  #               cov = t(round(exp2_res$CFA$cond1$ci_cov[1:3, 1:2], 0)))
  # # Tables
  # 
  # # Structure of resutls maethod x parameter
  # pt2 <- structure(temp, dim = c(2, ncol(temp)*2))
  # 
  # # Give names
  # pt3 <- rbind(c('bias', 'cov'), pt2) # column names
  # rownames(pt3) <- c(" ", # Empty space
  #                    "DURR", "IURR")

  # Generate base table
  pt4 <- xtable(pt3,
                align = c("r", rep("c", ncol(pt3))))
  
  # Modify Columns
  addtorow <- list()
  addtorow$pos <- list(0, 1)
  addtorow$command <- c(paste0(
    paste0('& \\multicolumn{2}{c}{', 
           c("f1", "f2", "f3"), 
           '}', 
           collapse=''), 
    '\\\\\n'),
    paste0(paste0('& \\multicolumn{6}{c}{', 
                  c("pm = .1, lv = 100"), 
                  '}', 
                  collapse=''), 
           '\\\\\n')
  )
  
  print(pt4, add.to.row = addtorow, include.colnames = F)
  print.xtable()
  
  Grade3 <- c("A","B","B","A","B","C","C","D","A","B",
              "C","C","C","D","B","B","D","C","C","D")
  Grade6 <- c("A","A","A","B","B","B","B","B","C","C",
              "A","C","C","C","D","D","D","D","D","D")
  Cohort <- table(Grade3, Grade6)
  Cohort
  xtable(Cohort)
  
  addtorow <- list()
  addtorow$pos <- list(0, 3)
  addtorow$command <- c("& \\multicolumn{4}{c}{Grade 6} \\\\\n",
                        "Grade 3 & A & B & C & D \\\\\n")
  print(xtable(Cohort), add.to.row = addtorow, include.colnames = FALSE)
  
  # Kable version
  
  # cond.seq <- which(1:8 %% 2 != 0) # odd conditions (p is the same)
  cond.seq <- c(1, 2)
  parm.seq <- 1:3 # which paramters do you want? The first 3?
  meth.seq <- 1:3 # which methods do you want?
  
  # Part 1: Join bias and CI for a given condition and analysis method
  
  exp2_res$conds
  
  cond.obj <- conds[cond.seq, ]
  
  pt1 <- lapply(cond.seq, function(x){
    cond.obj <- conds[x, ]
    target.obj <- exp2_res$CFA[[x]]
    temp <- rbind(Bias = t(target.obj$bias_per[parm.seq, 
                                               meth.seq+1]),
                  cov = t(round(target.obj$ci_cov[parm.seq, 
                                                  meth.seq], 0)))
    meth.names <- unique(rownames(temp))
    temp <- structure(temp, dim = c(max(meth.seq), length(meth.names)*2))
    colnames(temp) <- rep(c("bias", "cov"), length(meth.names))

    pt2 <- data.frame(cond = rep(paste0("pm = ", cond.obj$pm), each = nrow(temp)),
                      method = meth.names,
                      temp)
  })
  
  pt3 <- do.call(rbind, pt1)
  
  kbl(pt3, booktabs = T, align = "c") %>%
    add_header_above(c(" " = 1, " " = 1, "z1" = 2, "z2" = 2, "z3" = 2)) %>%
    pack_rows("lv = 10", 1, 3) %>%
    pack_rows("lv = 100", 4, 6) %>%
    collapse_rows(columns = 1)
  
  
  kbl(pt3, booktabs = T, align = "c", format = "latex") %>%
    add_header_above(c(" " = 1, " " = 1, "z1" = 2, "z2" = 2, "z3" = 2)) %>%
    pack_rows("lv = 10", 1, 3, latex_align = "c") %>%
    pack_rows("lv = 100", 4, 6, latex_align = "c") %>%
    collapse_rows(columns = 1)
  
  #############
  library(tables)
  mydf <- data.frame(rowFactor1 = sample(letters[1:2], 100, replace = TRUE), 
                     colFactor1 = sample(LETTERS[1:2], 100, replace = TRUE), 
                     x = rnorm(100), 
                     rowFactor2 = sample(1:2, 100, replace = TRUE), 
                     colFactor2 = sample(1:2, 100, replace = TRUE))
  
  tab1 <- tabular(Heading()*RowFactor(rowFactor2, spacing = 1, 
                                      levelnames = c("rowLabel1", "rowLabel2"))*
                    Heading()*RowFactor(rowFactor1, 
                                        levelnames = c("b1", "b2")) ~ 
                    Heading()*Factor(colFactor2, 
                                     levelnames = c("colLabel1", "colLabel2") )*
                    Heading()*Factor(colFactor1, 
                                     levelnames = c("a1", "a2"))*
                    Heading()*(x)*Heading()*(mean), 
                  data = mydf)
  
  Hmisc::latex(tab1)
  
#> lm scored data ####
  
  # # Bias
  # # High Factor loadings
  # lapply(exp2_res$lm,
  #        function(x) x$bias_per)[cindex_lh]
  # # Low factor loadings
  # lapply(exp2_res$lm,
  #        function(x) x$bias_per)[-cindex_lh]
  # 
  # # CI
  # lapply(exp2_res$lm,
  #        function(x) x$ci_cov)[cindex_lh]
  # lapply(exp2_res$lm,
  #        function(x) x$ci_cov)[-cindex_lh]

# # Summary Table -----------------------------------------------------------
# # This is a selected paramters version for summary paper presentation.
# 
#   col_id <- colnames(sum_exp1_sem$cond4$bias_raw)[c(12, 10, 2:9)]
#   indx   <- rownames(sum_exp1_sem$cond4$bias_raw)[c(1, 4, 
#                                                     7, 10, 
#                                                     13, 15, 25)]
#   sum_exp1_sem$cond4$cond
#   sum_exp1_sem$cond4$ci_cov[indx, ]
#   
#   sum_exp1_sem$cond4$cond
#   sum_exp1_sem$cond4$bias_per[indx, ]
#   
#   genTableEAM <- function(x){
#     # Generates section of the table for EAM turn in paper
#     store <- NULL
#     for (i in 1:length(indx)) {
#       store <- cbind(store,
#                      t(sum_exp1_sem[[x]]$bias_per[indx, col_id])[, i],
#                      t(sum_exp1_sem[[x]]$ci_cov[indx, col_id])[, i])
#     }
#     return(store)
#   }
#   
#   list_tables <- lapply(list(cond1=1,
#                              cond2=2,
#                              cond3=3,
#                              cond4=4), 
#                         genTableEAM)
#   
#   mt_table <- do.call(rbind, lapply(list_tables, rbind, NA))
# 
#   write.csv(mt_table, paste0("../output/", filename, "_table.csv") )
