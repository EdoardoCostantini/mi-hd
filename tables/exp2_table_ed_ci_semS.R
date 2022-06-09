### Title:    Analysis of results from experiment 2: euclidean distance for CIs
### Project:  Imputing High Dimensional Data
### Author:   Anonymized for peer review
### Created:  2020-07-09
### Notes:    wd should be the code folder

source("../tables/exp2_table_init.R")
exp2_res # check presence

# Table content -----------------------------------------------------------

  table.ed.obj <- round(t(sapply(exp2_res$semS,
                                 res_ed_ci)), 0)
  colnames(table.ed.obj) <- gsub(pattern = "\\_la", replacement = "", 
                                 colnames(table.ed.obj))
  colnames(table.ed.obj) <- gsub(pattern = "\\_", replacement = " ", 
                                 colnames(table.ed.obj))

# Format table ------------------------------------------------------------

  table.ed <- xtable(table.ed.obj, 
                     align = c("l",# 2nd slot is for actual rownames
                               rep("c", ncol(table.ed.obj)) ),
                     digits = 0,
                     caption = "Euclidean distance of the vector of 95\\% CIs"
  )

  addtorow         <- list()
  addtorow$pos     <- list(0, 0, 4, 4) # Add Condition 2 break
  addtorow$command <- c(# Add line break
    paste0('\\hline \n'),
    # factor loadings separation
    paste0(paste0('& \\multicolumn{',
                  ncol(table.ed)-1,
                  '}{c}{high $\\lambda$}', 
                  collapse=''), 
           '\\\\\n'),
    # Add line break
    paste0('\\hline \n'),
    # Condition Break
    paste0(paste0('& \\multicolumn{',
                  ncol(table.ed)-1,
                  '}{c}{low $\\lambda$}', 
                  collapse=''), 
           '\\\\\n')
  )
  
# Obtain latex code -------------------------------------------------------
  
  print(table.ed, 
        hline.after = c(0, 4,
                        nrow(table.ed.obj)), # repeating horizonta
        add.to.row = addtorow,
        include.rownames = TRUE,
        scalebox = 0.7)
