### Title:    Reading EVS data
### Author:   Edoardo Costantini
### Created:  2020-06-22

rm(list=ls())
library(foreign) # to import .dta data
library(labelled) # to extract valriables labels

# Read Data ---------------------------------------------------------------

  dat <- 1
  file_loc <- c("/Users/Work/Data/ZA7500_v3-0-0.dta/",
                "/Users/Work/Data/ZA7502_v1-0-0.dta/")[dat]
  
  file_name <- c("ZA7500_v3-0-0.dta",
                 "ZA7502_v1-0-0.dta")[dat]
  
  file_miss <- c("ZA7500_v3-0-0_missing.txt",
                 "ZA7502_v1-0-0_missing.txt")[dat]
  
  EVS2017 <- haven::read_dta(paste0(file_loc, file_name))
  dim(EVS2017)
  
  var_exclude <- c("v277", "mm_v277_fu")
  EVS2017 <- EVS2017[, -c(which( names(EVS2017) %in% var_exclude ))]
  EVS2017_GE_prepro <- as.data.frame(EVS2017[EVS2017$country == "276", ])
  EVS2017_GE <- as.data.frame(EVS2017[EVS2017$country == "276", ])
  
  names(EVS2017_GE)
  dim(EVS2017_GE)

# Read Missing data info --------------------------------------------------
  
  # Studying NA resons Replace NAs
  # Example of report missing types
  labelled::val_labels(EVS2017_GE[, "v171"])
  
  # Find all report missing type
  resp_cat <- lapply(EVS2017_GE, val_labels)
  mis_cat <- lapply(resp_cat, function(x) x < 0)
  unique(lapply(1:length(resp_cat), 
                function(x) resp_cat[[x]][mis_cat[[x]]]))
  
  # Read categories missing
  
  mis_codes <- read.table(paste0(file_loc, file_miss), sep	= "(")
  dim(mis_codes)
  
  mis_vars <- row.names(mis_codes)
    row.names(mis_codes) <- NULL
  mis_codes <- as.character(mis_codes[[1]])

  # Replace w/ NAs
  # Fix varnames
  mis_vars <- gsub(pattern = " ", 
                    replacement = "",
                    mis_vars)
  mis_codes <- mis_codes[-c(which( mis_vars %in% c("v277", "mm_v277_fu") ))]
  mis_vars <- mis_vars[-c(which( mis_vars %in% c("v277", "mm_v277_fu") ))]
  
  # get rid of period
  mis_codes <- gsub(pattern = "\\.", 
                    replacement = "",
                    mis_codes)
  
  # get rid of closing bracket
  mis_codes <- gsub(pattern = ")", 
                    replacement = "",
                    mis_codes)
  
  # Missings LOWEST THRU
  mis_codes <- gsub(pattern = "LOWEST THRU ", 
                    replacement = "",
                    mis_codes)
  
  # Get rid of negative values
  for (i in 1:10) {
    mis_codes <- gsub(pattern = as.character(-c(1:10))[i], 
                      replacement = "",
                      mis_codes)
  }
  
  mis_codes <- gsub(pattern = ", ", 
                    replacement = "",
                    mis_codes)
  mis_codes <- as.numeric(mis_codes)
  
  # Which variables in the original dataset should be worked
  mis_var_indx <- names(EVS2017_GE) %in% mis_vars

  # Replace w/ NAs
  for (i in 1:sum(mis_var_indx)) {
    # EVS2017_GE[, mis_var_indx][i] 
    # table(EVS2017_GE[, mis_var_indx][i])
    # Negative values to NAs
    # Other codes to NAs
    if(is.na(mis_codes[i])){
      NAindx <- EVS2017_GE[, mis_var_indx][i] < 0
      sum(EVS2017_GE[, mis_var_indx][i] < 0)
    } else {
      NAindx <- EVS2017_GE[, mis_var_indx][i] < 0 | 
        EVS2017_GE[, mis_var_indx][i] == mis_codes[i]
    }
    if(sum(NAindx) != 0){
      EVS2017_GE[NAindx, mis_var_indx][i] <- NA
    }
  }

# Study missingness ---------------------------------------------------------
  N <- nrow(EVS2017_GE)
  missPat <- mice::md.pattern(EVS2017_GE[, mis_var_indx], plot = FALSE)
  
## variablewise missing counts/proportions:
  missPat[nrow(missPat), ]
  missPat[nrow(missPat), -ncol(missPat)]
  
  percent_m  <- round(missPat[nrow(missPat), -ncol(missPat)]/N, 3)
  var_ll <- cbind(sapply(EVS2017_GE[, names(percent_m)], var_label))
  
  data.frame(
    percent_m = percent_m,
    label = var_ll
  )
  
  # Which vairables have large pm
  # 100% missings are country specific vairables not measures in Italy
  cbind(
    names(percent_m[percent_m == 1]),
    var_ll[percent_m == 1]
  )
  
  # pm > .8 & != 1
  cbind(
    percent_m[percent_m > .8 & percent_m != 1],
    var_ll[percent_m > .8 & percent_m != 1]
  )
  table(EVS2017_GE_prepro$v231b_r)
    # again due to non aplicability: if father was bonr in same country
    # this question is -3, not applicable
  table(EVS2017_GE_prepro$v248a)
  
  # pm > .35 & < .8
  cbind(
    percent_m[percent_m > .35 & percent_m < .8],
    var_ll[percent_m > .35 & percent_m < .8]
  )
  table(EVS2017_GE_prepro$v237) # gated question
  table(EVS2017_GE_prepro$v53) # gated question
  table(EVS2017_GE_prepro$v252_edulvlb) # gated question
  table(EVS2017_GE_prepro$v248) # gated question
  table(EVS2017_GE_prepro$v249) # gated question
  
  # pm < .35 & != 0
  cbind(
    percent_m[percent_m < .35 & percent_m != 0],
    var_ll[percent_m < .35 & percent_m != 0]
  )
  
  # 30% pm
  table(EVS2017_GE_prepro$v249) # gated question
    var_label(EVS2017_GE_prepro$v249)
  
  # 20% pm
  table(EVS2017_GE_prepro$v246_ISCO_2) # gated question
    var_label(EVS2017_GE_prepro$v246_ISCO_2)
  table(EVS2017_GE_prepro$v102) # not gated question
    var_label(EVS2017_GE_prepro$v102)
  
  # 10% pm
  table(EVS2017_GE_prepro$v31) # missing is due to nonresponse
    var_label(EVS2017_GE_prepro$v31)
  table(EVS2017_GE_prepro$v121) # missing is due to nonresponse
    var_label(EVS2017_GE_prepro$v121)
  
  # 5% 
  table(EVS2017_GE_prepro$v38) # missing is due to nonresponse
    var_label(EVS2017_GE_prepro$v38)
  
## Extract the patternwise missing:
  # How many missing variables per pattern
  missPat[-nrow(missPat), ncol(missPat)] # variables w/ missing in pattern
  
  length(missPat[-nrow(missPat), ncol(missPat)]) # total missing data patterns
  hist(sort(missPat[-nrow(missPat), ncol(missPat)]/ncol(missPat)),
       ylab = "")
  
## Compute covariance coverage:
  cc <- mice::md.pairs(EVS2017_GE[, mis_var_indx])$rr / nrow(EVS2017_GE[, mis_var_indx])
  unique(cc[cc < .80 & cc != 0])
  
### Study matrix desing
  table(EVS2017_GE_prepro$mm_matrix_group)
  table(EVS2017$mm_matrix_group)
  table(EVS2017$v13)
  table(EVS2017_GE_prepro$v41)
  

# Study distirbutions -----------------------------------------------------
  # Democratic likert-10
  mean(EVS2017_GE$v133_11c, na.rm = TRUE)
  dem_0 <- which(colnames(EVS2017_GE) == "v133")
  dem_1 <- which(colnames(EVS2017_GE) == "v144")
  dem_no <- which(grepl("_", colnames(EVS2017_GE)[dem_0:dem_1]))
  dem_it <- EVS2017_GE[, c(dem_0:dem_1)[-dem_no]]
  colnames(dem_it) <- paste0("dem", 1:ncol(dem_it))
  
  # Morals likert-10
  mor_0 <- which(colnames(EVS2017_GE) == "v149")
  mor_1 <- which(colnames(EVS2017_GE) == "v163")
  mor_it <- EVS2017_GE[, mor_0:mor_1]
  colnames(mor_it) <- paste0("mor", 1:ncol(mor_it))
  
  # Join data
  likert10 <- as.data.frame(cbind(mor_it, dem_it))
  
  # Means and variances
  mv <- data.frame(
    Means = apply(likert10, 2, mean, na.rm = TRUE),
    Vars = sapply(likert10, var, use = "pairwise.complete.obs")
  )
  round(mv, 3)
  colMeans(mv)

  # Range of covariances
  likert10_cov <- cov(likert10, use = "pairwise.complete.obs")
  mean(likert10_cov[lower.tri(likert10_cov)])
  range(likert10_cov[lower.tri(likert10_cov)])
  hist(likert10_cov[lower.tri(likert10_cov)])
  
  lapply(list(
    data.frame(
      Means = colMeans(likert10, na.rm = TRUE),
      Vars = sapply(likert10, var, use = "pairwise.complete.obs")
    ),
    cov(likert10, use = "pairwise.complete.obs"),
    cor(likert10, use = "pairwise.complete.obs")
  ), round, 3)
  
# Study CFA models --------------------------------------------------------

  FV <- c("v82", "v83", "v84")
  GR <- c("v72_DE", "v73_DE", "v74_DE", "v75_DE", 
          "v76_DE", "v77_DE", "v78_DE", "v79_DE")
  
  model <- '
    # Measurement Model
    GR =~ v72_DE + v73_DE + v74_DE + v75_DE + 
          v76_DE + v77_DE + v78_DE + v79_DE
    FV =~ v82 + v83 + v84
    
    # Latent relationships
    GR ~~ FV
  '
  
  fit <- cfa(model, data = EVS2017_GE, 
             std.lv = TRUE,
             missing = "fiml")
  
  # Check model fit
  summary(fit, fit.measures = TRUE, standardized = TRUE)
  parameterEstimates(fit, standardized = TRUE)
  

# Study Missingness of values scales --------------------------------------
  N <- nrow(EVS2017_GE)
  missPat <- mice::md.pattern(EVS2017_GE[, c(FV, GR)], plot = FALSE)
  
  missPat[nrow(missPat), ]
  missPat[95:nrow(missPat), -ncol(missPat)]
  
  pm  <- 1-round(missPat[nrow(missPat), -ncol(missPat)]/N, 3)
  var_ll <- cbind(sapply(EVS2017_GE[, names(pm)], var_label))
  
  data.frame(
    percent_m = pm,
    label = var_ll
  )
  
  
  
  
  