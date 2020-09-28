### Title:    Analysis of results from experiment 2: factor loadings table
### Project:  Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-07-09
### Notes:    wd should be the code folder

  source("../tables/exp2_table_init.R")
  exp2_res # check presence
  
# Table content -----------------------------------------------------------

  # Same pm, different lv
  cond.seq <- 3:4 # condition of interest (rows of exp2_res$conds)
  parm.seq <- 1:5 # parameters of interest (items of first latent variable)
  mi.meth.tab <- 1:8 # mi methods
  si.meth.tab <- 9:11 # not mi methods
  conds <- exp2_res$conds
  
  # Structure results for MI methods for selected conditions, parameters, and methods
  mi_bici <- lapply(cond.seq, function(x){
    # x <- cond.seq[1]
    # Extract Values of interest
    cond.obj   <- conds[x, ]
    target.obj <- exp2_res$CFA[[x]]
    
    # Shape display
    mi.bici <- rbind(Bias = t(target.obj$bias_per[parm.seq, 
                                                  mi.meth.tab+1]),
                     cov = t(round(target.obj$ci_cov[parm.seq, 
                                                     mi.meth.tab], 0)))
    meth.names <- gsub(pattern = "\\_la", replacement = "", 
                       unique(rownames(mi.bici)))
    meth.names <- gsub(pattern = "\\_", replacement = " ", 
                       meth.names)
    
    mi.bici <- structure(mi.bici, dim = c(length(mi.meth.tab), length(parm.seq)*2))
    
    # Produce output
    mi_bici.out <- cbind(meths = meth.names, mi.bici)
    return(mi_bici.out)
  })
  
  # Join results
  mi_bici <- do.call(rbind, mi_bici)
  
  # Structure results for SI methods for selected conditions, parameters, and methods
  
  si_bici <- lapply(cond.seq[1], # only for 1 for pm
                    function(x){
    
    # Extract Values of interest
    cond.obj   <- conds[x, ]
    target.obj <- exp2_res$CFA[[x]]
    temp <- rbind(Bias = t(target.obj$bias_per[parm.seq, 
                                               si.meth.tab+1]),
                  cov = t(round(target.obj$ci_cov[parm.seq, 
                                                  si.meth.tab], 0)))
    meth.names <- unique(rownames(temp))
    temp <- structure(temp, dim = c(length(si.meth.tab), length(parm.seq)*2))
    si_bici.out <- cbind(meths = meth.names, temp)
    return(si_bici.out)
  })
  
  si_bici <- do.call(rbind, si_bici)
  
  content <- rbind(si_bici, mi_bici)
  colnames(content) <- c("methods", rep(c("Bias", "CR"), length(parm.seq)))
  
# Format table ------------------------------------------------------------

  table.fl <- xtable(x = content, 
                     caption = "Factor loadings of the first 5 items for the first latent variable",
                     align = c("", # 1st slot is for R object rownames (not displayed)
                               "l",# 2nd slot is for actual rownames
                               rep(c("|c", "c"), (ncol(mi_bici)-1)/2 ))
  )
  
  # Fix Column names and the like
  addtorow         <- list()
  addtorow$pos     <- list(-1, # Add latent factors header
                           -1, # line break
                           3,  # line break
                           3, # Add Condition 1 break
                           11, # line break
                           11) # Add Condition 2 break
  addtorow$command <- c(
    # Factor loading name
    paste0(paste0('& \\multicolumn{2}{c}{', 
                  paste0("$\\lambda_", parm.seq, "$"), 
                  '}', 
                  collapse=''), '\\\\\n'),
    # Add line break
    paste0('\\hline \n'),
    # Add line break
    paste0('\\hline \n'),
    # Condition Break
    paste0(paste0('\\multicolumn{',
                  ncol(table.fl),
                  '}{c}{', 
                  paste0("pm = ", conds[cond.seq[1], ]$pm,
                         ", lv = ", conds[cond.seq[1], ]$lv),
                  '}', 
                  collapse=''), 
           '\\\\\n'),
    # Add line break
    paste0('\\hline \n'),
    # Condition Break 
    paste0(paste0('\\multicolumn{',
                  ncol(table.fl),
                  '}{c}{', 
                  paste0("pm = ", conds[cond.seq[2], ]$pm,
                         ", lv = ", conds[cond.seq[2], ]$lv),
                  '}', 
                  collapse=''), 
           '\\\\\n')
  )
  
# Obtain latex code -------------------------------------------------------

  print(table.fl, 
        hline.after = c(0, 
                        3,
                        11,
                        nrow(table.fl)), # repeating horizonta
        add.to.row = addtorow, 
        include.colnames = TRUE,
        include.rownames = FALSE,
        scalebox = 0.75)
  