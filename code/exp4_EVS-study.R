### Title:    Reading EVS data for resampoling experiment and studying how 
###           it is structured
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

# Study Integrated --------------------------------------------------------

  # Variables Types
  # Structural missings: religiosity
  int.df[int.df$v51 == 2, "v53"]
  int.df[int.df$v51 == 1, "v53"] 
  # these are missings because question v53 is asked only if v51 != 1
  
  # Country specific questions
  colnames(int.df)
  
  int.df[, c("v24", "v24a_IT", "v24b_IT")]
  
  # Italy experimental formualtion of v24
  int.df[country_id == 380, c("v24", "v24a_IT", "v24b_IT")]
  
  # # Use v24a and v24b to substitute NA in v24 if possible
  # for (i in 1:nrow(int.df)) {
  #   if(is.na(int.df[i, "v24"])){
  #     int.df[i, "v24"] <- int.df[i, "v24a_IT"]
  #     if(is.na(int.df[i, "v24"])){
  #       int.df[i, "v24"] <- int.df[i, "v24b_IT"]
  #     }
  #   }
  # }
  
  tab_country(int.df, int.df$country, "v224")
  tab_country(int.df, int.df$country, "v224_DK")

# What to keep what to throw ----------------------------------------------

  id  <- "id_cocas"
  ord <- paste0("v", c(1:8, 32:39, 46:50, 63:70, 72:84,
                       97:107, 115:168, 170:172,
                       176:203, 205:224, 226,
                       240, 242, 247, 267:274, 280, 
                       c("174_LR",  "239a", "239b", "261_ppp")))
  nom <- paste0("v", 52)
  
  # Income
  tab_country(int.df, int.df$country, "v261")
  int.df$v261_ppp  # income corrected for purchasing power parity (PPP)
  inc <- "v261_ppp"
  
  # Age
  int.df$v226 
  int.df$age  # at time of questionnaire
  int.df$age_r
  int.df$age_r2
  int.df$age_r3
  
  age <- "age"
  
  # Education
  # in the codebook table with all names related to education
  int.df$v243_ISCED_1 # respondant
  int.df$v252_ISCED_1 # partner (gated!)
  int.df$v262_ISCED_1 # father
  int.df$v263_ISCED_1 # mother
  
  edu <- c("v243_ISCED_1", 
           # "v252_ISCED_1", # Gated
           "v262_ISCED_1", 
           "v263_ISCED_1")
  
  # Job Profession
  # in the codebook table with all names related to job
  
  # Income
  int.df$v261_ppp
  
  # Recoded variables (e.g. nuymber of children)
  head(int.df[, c("v239a", "v239b", "v239_r")], 10)
  
  # Political Parties
  tab_country(int.df, int.df$country, "v174_LR")
  tab_country(int.df, int.df$country, "v175_LR")
  as.numeric(int.df$v174_LR)
  sort(unique(as.numeric(int.df$v174_LR)))
  int.df$v175_LR
  
  # Structured missinges
  str_mis <- c("v229", "v53", "v231b", "v231b_r", "v233b", "v233b_r", 
               "v235", "v236", "v237", "v241", "v248a", "v249", "v250", 
               "v251b", "v251b_r")
  
  # Country of birth vairables
  
  # Keep selected variables (fl = filtered)
  int.dt.fl <- int.df[, c(id, ord, edu, inc, age)]
  
  str_mis[which(str_mis %in% colnames(int.dt.fl))] # no more left

# Reorganize --------------------------------------------------------------
  list.df <- lapply(list(int.df = int.df, mad.df = mad.df), clean_up,
                    id  = "id_cocas",
                    country = "country",
                    age = "age",
                    inc = "v261_ppp",
                    ord = paste0("v", c(1:8, 32:39, 46:50, 63:70, 72:84,
                                        97:107, 115:168, 170:172,
                                        176:203, 205:224, 226,
                                        240, 242, 247, 267:274, 280, 
                                        c("174_LR",  "239a", "239b", "261_ppp"))),
                    edu = c("v243_ISCED_1", 
                            # "v252_ISCED_1", # Gated
                            "v262_ISCED_1", 
                            "v263_ISCED_1"))
  
  sapply(list.df, function(x) table(x$country))
  
  x <- list.df$int.df
  
  length(which(list.df$int.df[, "id"] %in% list.df$mad.df[, "id"])) # duplicate cases
  sum(int.df$fduplicate)
  
  list.df$int.df$id

# Missingness -------------------------------------------------------------

  # Transform negative values in NA (R value)
  int.dt.fl[int.dt.fl < 0] <- NA
  
  # Check how many fully observed variables you have
  mdPat_small <- md.pattern(int.dt.fl[, 1:4], plot = FALSE)
  c(nMisVar = mdPat_small[1, ncol(mdPat_small)],
    nObsCas = as.numeric(rownames(mdPat_small)[1]))
  
  sum((rowSums(is.na(int.dt.fl)) == 0))
  mdPat <- md.pattern(int.dt.fl, plot = FALSE)
  c(nMisVar = mdPat[1, ncol(mdPat)],
    nObsCas = as.numeric(rownames(mdPat)[1]))

# Study missingness ---------------------------------------------------------
  N <- nrow(int.dt.fl)
  missPat <- mice::md.pattern(int.dt.fl, plot = FALSE)
  
  ## variablewise missing counts/proportions:
  missPat[nrow(missPat), ]
  missPat[nrow(missPat), -ncol(missPat)]
  
  percent_m  <- round(missPat[nrow(missPat), -ncol(missPat)]/N, 3)
  var_ll <- cbind(sapply(int.dt.fl[, names(percent_m)], var_label))
  
  data.frame(
    percent_m = percent_m,
    label = var_ll
  )
  
  # Which vairables have large pm
  
  # pm > .8 & != 1
  cbind(
    percent_m[percent_m > .8 & percent_m != 1],
    var_ll[percent_m > .8 & percent_m != 1]
  )
  
  # pm > .35 & < .8
  cbind(
    percent_m[percent_m > .35 & percent_m < .8],
    var_ll[percent_m > .35 & percent_m < .8]
  )
  
  # pm < .35 & != 0
  cbind(
    percent_m[percent_m < .35 & percent_m != 0],
    var_ll[percent_m < .35 & percent_m != 0]
  )

## Extract the patternwise missing:
# How many missing variables per pattern
missPat[-nrow(missPat), ncol(missPat)] # variables w/ missing in pattern
length(missPat[-nrow(missPat), ncol(missPat)]) # total missing data patterns
hist(sort(missPat[-nrow(missPat), ncol(missPat)]/ncol(missPat)),
     ylab = "")

## Compute covariance coverage:
cc <- mice::md.pairs(int.dt.fl)$rr / nrow(int.dt.fl)
unique(cc[cc < .80 & cc != 0])


# Matrix Design -----------------------------------------------------------

rm(list=ls())
source("./init_general.R")
source("./exp4_init.R")

# Load Data Variables
  data_source <- readRDS("../data/exp4_EVS2017_full.rds")$full
  # Clean Variable Names
  EVS_vars <- gsub("_(.*)", "", colnames(data_source))

# Load Variable-Block legend
  X <- read.csv("~/Data/EVS2017/documents/ZA7500_vars_blcok_membership.csv",
           header = TRUE)
  class(X) # check type
  
  # Clean Block names
  X$Variable <- tolower(trimws(as.character(X$Variable)))
  
# Obtain Variables Block Classification
  var_class <- data.frame(Variable = X$Variable[X$Variable %in% EVS_vars],
                          Block = X$Block[X$Variable %in% EVS_vars])
  extra <- data.frame(Variable = EVS_vars[!EVS_vars %in% X$Variable],
                      Block = c("Core","Core","Core","Core","Core"))
  var_class_f <- rbind(var_class, extra)

# Clean Data
  var_class_f$Variable <- as.character(var_class_f$Variable)
  var_class_f$Block <- droplevels(var_class_f$Block)
  var_class_f$Block <- factor(var_class_f$Block, levels = c("Core", "A", "B", "C", "D"))
  
  table(var_class_f$Block)
  
## Data ------------------------------------------------------------------ ##
# Gen one fully-obs data
Xy <- data_source[sample(1:nrow(data_source),
                         cond$n,
                         replace = TRUE), ]
