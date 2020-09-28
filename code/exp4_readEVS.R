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
  int.df$country_org <- int.df$country #original country variable
  int.df$country <- factor(int.df[, "country"],
                           levels = val_labels(int.df[, "country"]),
                           labels = names(val_labels(int.df[, "country"])))

# Decisions ---------------------------------------------------------------
  
  # Subset data
  subsample <- c(1, 4)
  countries  <- c("France", "Germany", "Italy", "Netherlands")
  
  int.df <- int.df %>%
    filter(mm_select_sample %in% subsample) %>%
    filter(country %in% countries)
  
  # Variables Selection
  id  <- "id_cocas"
  ord <- paste0("v", c(1:8, 32:39, 46:50, 53:54, 63:70, 72:84,
                       97:107, 115:168, 170:172, 176:203, 205:224, 
                       240, 242, 267:274, 280, 
                       c("174_LR", "175_LR", "239_r", "239a", "239b")))
  dic <- paste0("v", c(9:31,40:45,51,57:61,71,85:95,
                       112,169,225,227,230,232,236,
                       248,259,260))
  nom <- paste0("v", c(52, 62, 56, # religiosity
                       108:111, 113:114, 234, 244, 238,
                       264, 265, 266, "276_r")) 
  
  # I will make sure the following variables are excluded
  excl <- c("v27", # not really a variable in the 
            # countries we are investigating (98% of responses 
            # was same category)
            "v239a", "v239b",
            # Gated
            "v53", "v247", "v248", "v252_ISCED_1", "v170",
            # Religiosity related
            "v51", "v52",
            # respondend related questions to be discarded after use
            "v228b_r",
            # father related question to discard after use
            "v230", "v231b", "v231b_r",
            # mother related question to discard after use
            "v232", "v233b", "v233b_r"
  )
  
# Prepare Variables -------------------------------------------------------
  
# Political Parties
  pol <- c("v175_LR", "v174_LR")
  
# Age
  age <- "age"
  
# Education
  edu <- c("v243_ISCED_1", "v262_ISCED_1", "v263_ISCED_1")
  
# Job / SES
  SES <- "v246_egp"

# Income
  inc <- "v261_ppp"
  
# Religiosity
  den <- c("v51", # do you belong
           "v52") # which denomination do you belong to?
  rel <- c("v54", "v55") # apart from ... how often services? present / future
  
  # Combine Information from the two variables
  v51v52_comb <- case_when(int.df[, den[1]] == 1 ~ int.df[, den[2]],
                        # If you belong to a denomination, to which?
                        int.df[, den[1]] == 2 ~ 0)
                        # If you do not, then have value 0
  val_labels(v51v52_comb) <- c("No Religion" = 0, 
                               val_labels(v51v52_comb)[val_labels(v51v52_comb) > 0])
  int.df$v51v52_comb <- v51v52_comb
  
# Marital Status
  mart <- "v234"
  
# Sex
  sex <- "v225"
  int.df$v225 <- recode(int.df$v225, "2" = 0)
  
# Union Membership
  um <- "v11"
  
# Countries of birth (repsondent, father, mother)
  # Map country names to codes
  ISO3116 <- data.frame(name = names(val_labels(int.df[, "v231b_r"])), 
                        code = as.numeric(val_labels(int.df[, "v231b_r"])))
  
  # Interviwee country of birth
  re_intw  <- ISO3116$name[match(int.df[, "country_org"], ISO3116$code)] # interview
  re_alter <- ISO3116$name[match(int.df[, "v228b_r"], ISO3116$code)] # birth
  re_cb <- case_when(re_alter == "not applicable" ~ re_intw, # country of birth
                     re_alter != "not applicable" ~ re_alter)
  
  
  # Father's country
    fa_alter <- ISO3116$name[match(int.df[, "v231b_r"], ISO3116$code)]
    fa_cb    <- case_when(int.df[, "v230"] == 1 ~ re_intw,
                          int.df[, "v230"] == 2 ~ fa_alter)
    # Add to dataset
    int.df$fa_cb <- fa_cb
    
    father_index <- c("fa_cb")
    
  # Mother's country
    ma_alter <- ISO3116$name[match(int.df[, "v233b_r"], ISO3116$code)]
                        # Respondent born in country of interview -> keep that value
    ma_cb <- case_when(int.df[, "v232"] == 1 ~ re_intw,
                        # REspondent born in other country -> check 
                       int.df[, "v232"] == 2 ~ ma_alter)
    # Add to dataset
    int.df$ma_cb <- ma_cb

# Define Vector of variable names
  # Prepare vector of column names to keep
  var_indx <- unique(c(id,
                       "country", 
                       gtools::mixedsort(c(ord, dic, nom, pol, age, mart, 
                                           sex, um, edu, SES, inc, 
                                           "v51v52_comb", rel)),
                       "ma_cb", "fa_cb"))
  
  var_vec <- var_indx[!var_indx %in% excl]
  
# Get desired data --------------------------------------------------------

  # Perform Variable Selection
  EVS2017 <- int.df %>%
    select(all_of( var_vec ))
  dim(EVS2017)
  # all_of gets rid of vector of names ambiguity
  
# Dealing with missing data -----------------------------------------------

  # Check for missing values differently coded and give them a 
  # known missing value code (-2)
  NA_map <- data.frame(
    codes = c(6, 66, 0, 7, 3, 9999),
    type = c("does not apply to me (spontaneous)", 
             "does not apply to me (spontaneous)",
             # I'm treating it as a value to impute based on other people 
             "no formal education", 
             # what do other similar people have? i'm going to impute
             "not allowed to vote",
             # what do other similar people have? i'm going to impute
             "other answer (code if volunteered only!)",
             "no answer")
  )

  for (v in 1:ncol(EVS2017)) {
    lab_index <- 
      val_labels(EVS2017[, v]) %in% NA_map$codes & 
      names(val_labels(EVS2017[, v])) %in% NA_map$type
    EVS2017[, v][EVS2017[, v] == val_labels(EVS2017[, v])[lab_index]] <- NA
  }
  
  # Assign NA to all coded as missings
  col_indx <- var_vec[!var_vec %in% c("fa_cb", "ma_cb", "country")]
  EVS2017[, col_indx][EVS2017[, col_indx] < 0] <- NA
  
# Transform Vairables to proper R objects ---------------------------------
# i.e. get rid of dang haven_labelled

  # Transfrom variables in class that I want
  EVS2017_cl <- EVS2017 %>%
    mutate(id_cocas = as.character(id_cocas)) %>%
    mutate_if(colnames(EVS2017) %in% id, as.character) %>%
    mutate_if(colnames(EVS2017) %in% "country", as.factor) %>%
    mutate_if(colnames(EVS2017) %in% age, as.numeric)  %>%
    mutate_if(colnames(EVS2017) %in% ord, as.numeric)  %>%
    mutate_if(colnames(EVS2017) %in% dic, as.factor)  %>%
    mutate_if(colnames(EVS2017) %in% nom, as.factor)   %>%
    mutate_if(colnames(EVS2017) %in% "v51v52_comb", as.factor)  %>%
    mutate_if(colnames(EVS2017) %in% rel, as.numeric)  %>%
    mutate_if(colnames(EVS2017) %in% edu, as.numeric)  %>%
    mutate_if(colnames(EVS2017) %in% inc, as.numeric)  %>%
    mutate_if(colnames(EVS2017) %in% pol, as.numeric)  %>%
    mutate_if(colnames(EVS2017) %in% SES, as.factor)
  
  # Fix some things
  EVS2017_cl$country <- droplevels(EVS2017_cl$country)
  levels(EVS2017_cl$v51v52_comb) <- names(val_labels(EVS2017$v51v52_comb))
  levels(EVS2017_cl$v246_egp) <- names(val_labels(EVS2017$v246_egp))[-c(1:7, 19)]
  
  cbind(
    table(EVS2017_cl$v246_egp),
    table(EVS2017$v246_egp)
  )

  # Make ID column as row name
  
  EVS2017_cl <- data.frame(EVS2017_cl, row.names = 1)[, -c(249, 250)]
  sapply(EVS2017_cl, class)
  
  # Imputation with PMM
  N <- nrow(EVS2017_cl)
  missPat <- mice::md.pattern(EVS2017_cl, plot = FALSE)
  
  ## variablewise missing counts/proportions:
  missPat[nrow(missPat), ]
  
  missPat[nrow(missPat), -ncol(missPat)]
  percent_m <- round(missPat[nrow(missPat), -ncol(missPat)]/N, 3)
  percent_m["v276_r"]
  predMat <- quickpred(EVS2017_cl, mincor = .3)
  
  # Check out which variables are used for which imputation models
  store <- data.frame(imp = rep(NA, ncol(predMat)),
                      pred = rep(NA, ncol(predMat)))
  
  for (v in 1:ncol(predMat)) {
    store[v, ] <- c(rownames(predMat)[v],
                    paste0(names(which(predMat[v, ] != 0)), collapse = ", "))
  }
  
  # # Perform Convergence check
  # mids_pmm <- mice::mice(EVS2017_cl, 
  #                        m = 5, maxit = 200,
  #                        predictorMatrix = predMat,
  #                        method = "pmm")
  # saveRDS(mids_pmm, "../output/exp4_ccheck_mids_pmm.rds")
  # 
  # # CHECK Imputations for EACH VAR TYPE
  # EVS2017_full <- complete(mids_pmm)
  # sapply(EVS2017_full, class)
  # 
  # lapply(var_list[5], function(x){
  #   list((EVS2017 %>%
  #            select(x))[1:20, 1],
  #         (EVS2017_full %>%
  #            select(x))[1:20, 1])
  # })
  # 
  # var_target <- names(which(missPat[nrow(missPat), ] != 0)) # target of miss
  # plot(mids_pmm, y = var_target[var_target %in% c(var_list$ord[21:30])])
  # 
  # # Perform Single imputation to generate desired data
  # mids_final <- mice::mice(EVS2017_cl, 
  #                          m = 1, maxit = 1e2,
  #                          predictorMatrix = predMat,
  #                          method = "pmm")
  
  # Fake version
  mids_final <- mice::mice(EVS2017_cl, 
                           m = 1, maxit = 1,
                           predictorMatrix = predMat,
                           method = "pmm")

  mids_final$loggedEvents
  
  EVS2017_psuedo <- complete(mids_final)
# 

  
  # Save Results of cleaning and imputation
  out <- list(orig = EVS2017_cl,
              full = EVS2017_psuedo,
              note = paste0("2017 EVS dataset for selected countries and ",
                            "variables, ready for resampling study.") 
              )
  
  saveRDS(out, "../data/exp4_EVS2017_full.rds")
  