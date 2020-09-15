### Title:    Reading EVS data for resampoling experiment
### Author:   Edoardo Costantini
### Created:  2020-06-22

rm(list=ls())
library(foreign) # to import .dta data
library(labelled) # to extract valriables labels
library(mice)
library(dplyr)

source("./functions_EVS.R")

# Read Data ---------------------------------------------------------------

  file_loc <- "/Users/Work/Data/EVS2017/data/"
  file_int <- "ZA7500_v3-0-0" # integrated data
  file_mat <- "ZA7502_v1-0-0" # matrix desing version
  miss_tag <- "_missing.txt"
  
  # Read Integrated data
  int.dt <- haven::read_dta(paste0(file_loc, file_int, ".dta"))
  dim(int.dt)
  int.df <- as.data.frame(int.dt)
  int.df$country <- factor(int.df[, "country"], 
                           levels = val_labels(int.df[, "country"]), 
                           labels = names(val_labels(int.df[, "country"])))
  
  mad.dt <- haven::read_dta(paste0(file_loc, file_mat, ".dta"))
  dim(mad.dt)
  mad.df <- as.data.frame(mad.dt)
  mad.df$country <- factor(mad.df[, "country"], 
                           levels = val_labels(mad.df[, "country"]), 
                           labels = names(val_labels(mad.df[, "country"])))

# Decisions ---------------------------------------------------------------
  
  # Choosing subsamples
  lapply(list(int.dt = int.dt$mm_select_sample,
              mad.dt = mad.dt$mm_select_sample), table)
  val_labels(int.dt$mm_select_sample)
  subsample <- 1
    # Conclusion: go w/ 1, regular interviewer-administred group
  
  # Countries
  val_labels(int.dt$country)
  lapply(list(int.dt = int.dt$country,
              mad.dt = mad.dt$country), table)
    # Conclusion: choose EU founders
  countries  <- c("Belgium", "France", "Germany", "Italy", "Luxembourg", "Netherlands")
  
  # Variables Selection
  id  <- "id_cocas"
  ord <- paste0("v", c(1:8, 32:39, 46:50, 63:70, 72:84,
                       97:107, 115:168, 170:172,
                       176:203, 205:224,
                       240, 242, 267:274, 280, 
                       c("174_LR", "239_r", "239a", "239b")))
  dic <- paste0("v", c(9:31,40:45,51,57:61,71,85:95,
                       112,169,225,227,230,232,248,259,260))
  nom <- paste0("v", c(52, 62,   # religiosity
                       108:111, 113:114, 234, 238)) 
  
  gated <- c("v52", "v53", 
             "v247", "v248", 
             "v252_ISCED_1",
             "v170")
    # I will exclude these vairables for sure as they would imply a complex imputation
    # set up that we are not really interested in including at the moment
  
  excl <- c("v27", # not really a variable in the 
                   # countries we are investigating (98% of responses 
                   # was same category)
            "v239a", "v239b")
  
  # Income
  inc <- "v261_ppp"
  
  # Age
  age <- "age"
  
  # Education
  edu <- c("v243_ISCED_1", 
           "v252_ISCED_1", # Gated
           "v262_ISCED_1", 
           "v263_ISCED_1")
  
  # Job/Profession/SES
  # Options:
  # int.dt$v246_ISEI
  # int.dt$v246_SIOPS
  SES <- "v246_ISEI"
    # Both numerical w/ respeact to others that are categorical
    # No idea on which is better
  
  # Income
  inc <- "v261_ppp"
  
  # Political Parties
  pol <- "v175_LR"
  
  # Prepare vector of column names to keep
  var_indx <- list(id, age, ord, dic, nom, edu, inc, pol, SES)
    names(var_indx) <- c("id", "age", "ord", "dic", "nom", "edu", "inc", "pol", "SES")
  var_list <- sapply(var_indx, function(x) x[!x %in% c(gated, excl)])
  var_vec  <- do.call(c, var_list)
    names(var_vec) <- NULL

# Get desired data --------------------------------------------------------

  # Perform Variable Selection
  EVS2017 <- int.df %>%
    filter(mm_select_sample == subsample) %>%
    filter(country %in% countries) %>%
    select(all_of( var_vec ))
      # all_of gets rid of vector of names ambiguity

  # Check for missing values differently coded and give them a 
  # known missing value code (-2)
  NA_codes <- c(6, 66, 0, 7, 3, 9999)
  NA_type <- c("does not apply to me", 
               "does not apply to me",
                # I'm treating it as a value to impute based on other people 
               "no formal education", 
                # what do other similar people have? i'm going to impute
               "not allowed to vote",
                # what do other similar people have? i'm going to impute
               "other answer",
               "no answer")
  
  for(v in seq_along(NA_codes)){
    outCount <- EVS2017 %>%
      summarise_each(funs(sum(. == NA_codes[v])))
    outNames <- names(outCount)[outCount > 0]
    outNames <- outNames[outNames %in% var_vec]
    outLabel <- sapply(EVS2017[, outNames], val_labels) # plausible NAs
    outMiss  <- names(outLabel[grep(NA_type[v], outLabel)])
    
    for (j in seq_along(outMiss)) {
      logi.indx <- EVS2017[, outMiss[j]] == NA_codes[v]
      EVS2017[logi.indx, outMiss[j]] <- -2 # assign "no an
    }
  }
  
  # Assign NA to all coded as missings
  EVS2017[EVS2017 < 0] <- NA
  
  # Transfrom variables in class that I want
  EVS2017_cl <- EVS2017 %>%
    mutate(id_cocas = as.character(id_cocas)) %>%
    mutate_if(colnames(EVS2017) %in% var_list$id, as.character) %>%
    mutate_if(colnames(EVS2017) %in% var_list$age, as.numeric)  %>%
    mutate_if(colnames(EVS2017) %in% var_list$ord, as.numeric)  %>%
    mutate_if(colnames(EVS2017) %in% var_list$dic, as.factor)  %>%
    mutate_if(colnames(EVS2017) %in% var_list$nom, as.factor)   %>%
    mutate_if(colnames(EVS2017) %in% var_list$edu, as.numeric)  %>%
    mutate_if(colnames(EVS2017) %in% var_list$inc, as.numeric)  %>%
    mutate_if(colnames(EVS2017) %in% var_list$pol, as.numeric)  %>%
    mutate_if(colnames(EVS2017) %in% var_list$SES, as.numeric)
  
  sapply(EVS2017_cl, class)
  
  # Imputation with PMM
  N <- nrow(EVS2017_cl)
  missPat <- mice::md.pattern(EVS2017_cl, plot = FALSE)
  
  ## variablewise missing counts/proportions:
  missPat[nrow(missPat), ]
  
  missPat[nrow(missPat), -ncol(missPat)]
  percent_m <- round(missPat[nrow(missPat), -ncol(missPat)]/N, 3)
  
  predMat <- quickpred(EVS2017_cl, mincor = .3)
  
  # Force all nominal variables that are still coded as such to be included
  predMat[, colnames(predMat) %in% nom] <- 1
  
  # Perform Convergence check
  mids_pmm <- mice::mice(EVS2017_cl, 
                         m = 5, maxit = 200,
                         predictorMatrix = predMat,
                         method = "pmm")
  saveRDS(mids_pmm, "../output/exp4_ccheck_mids_pmm.rds")
  
  # CHECK Imputations for EACH VAR TYPE
  EVS2017_full <- complete(mids_pmm)
  sapply(EVS2017_full, class)
  
  lapply(var_list[5], function(x){
    list((EVS2017 %>%
             select(x))[1:20, 1],
          (EVS2017_full %>%
             select(x))[1:20, 1])
  })
  
  var_target <- names(which(missPat[nrow(missPat), ] != 0)) # target of miss
  plot(mids_pmm, y = var_target[var_target %in% c(var_list$ord[21:30])])
  
  # Perform Single imputation to generate desired data
  mids_final <- mice::mice(EVS2017_cl, 
                           m = 1, maxit = 1e2,
                           predictorMatrix = predMat,
                           method = "pmm")
  EVS2017_psuedo <- complete(mids_pmm)
  
  out <- list(dt = EVS2017_psuedo,
              note = paste0("2017 EVS dataset for selected countries and ",
                            "variables, ready for resampling study.") 
              )
  
  saveRDS(out, "../data/exp4_EVS2017_full.rds")
  