### Title:    Reading EVS data for resampoling experiment
### Author:   Edoardo Costantini
### Created:  2020-06-22
  
  rm(list=ls())
  library(foreign) # to import .dta data
  library(labelled) # to extract valriables labels
  library(mice)
  library(dplyr)
  library(forcats) # for fct_collapse() function to recode factors
  
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
  # int.df$country <- factor(int.df[, "country"],
  #                          levels = val_labels(int.df[, "country"]),
  #                          labels = names(val_labels(int.df[, "country"])))
  length(val_labels(int.df$country))

# Decisions ---------------------------------------------------------------
  
  # Subset data
  subsample <- c(1, 4)
  # countries  <- c("France", "Germany", "Italy", "Netherlands")
  countries  <- c(France = 250, Germany = 276, Italy = 380, Netherlands = 528)
  
  int.df <- int.df %>%
    filter(mm_select_sample %in% subsample) %>%
    filter(country %in% countries)
  colnames(int.df)
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
  # int.df$v225 <- recode(int.df$v225, "2" = 0)
  
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
  EVS2017 <- dplyr::select(int.df, var_vec)
  dim(EVS2017)
  
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
  
  # Perform Convergence check
  mids_pmm <- mice::mice(EVS2017_cl,
                         m = 5, maxit = 200,
                         predictorMatrix = predMat,
                         method = "pmm")
  saveRDS(mids_pmm, "../convergence/exp4_ccheck_mids_pmm.rds")

  var_target <- names(which(missPat[nrow(missPat), ] != 0))
  plot(mids_pmm, y = var_target[var_target %in% var_vec[1:20]])
  
  # Final Version
  # mids_final <- mice::mice(EVS2017_cl, 
  #                          m = 1, maxit = 1e2, # NEED TO RUN SERIOUSLY
  #                          predictorMatrix = predMat,
  #                          method = "pmm")
  # mids_final$loggedEvents
  # EVS2017_psuedo <- complete(mids_final)
  
  # or simply save 1 dataset from the 5 you have in the mids_pmm output
  mids_pmm <- readRDS("../convergence/exp4_ccheck_mids_pmm.rds")
  EVS2017_psuedo <- complete(mids_pmm, action = 1)
  
## Fix factors ------------------------------------------------------------- ##
  
  EVS2017_final <- sapply(colnames(EVS2017_psuedo), 
                          fix.factor, 
                          dt_h = EVS2017,
                          dt_imp = EVS2017_psuedo,
                          simplify = FALSE,
                          USE.NAMES = TRUE)
  EVS2017_final <- do.call(data.frame, EVS2017_final)
  
## Fix Low variance variables ---------------------------------------------- ##
  # EVS2017_cl <- readRDS("../data/exp4_EVS2017_full.rds")$orig
  # EVS2017_final <- readRDS("../data/exp4_EVS2017_full.rds")$full

  # Generate model matrix
  mm <- model.matrix(~., EVS2017_final)[, -1]

  # Get rid of constants (might happen for 1 unpopular dummy code)
  const <- names(which(apply(mm, 2, var) == 0))
  mm  <- mm[, !colnames(mm) %in% const]
  
  # Find Dummies that have 95% of observations in 1 category
  tabular <- apply(mm, 2, table)
  
  # Select only dummies
  tabular <- tabular[sapply(tabular, length) == 2]
  
  # Vector of dummy names to discard / things to collapse
  dum.disc <- lapply(tabular, function(x) {
    x[1] / sum(x) > .95 | x[1] / sum(x) < 1-.95
  })
  dum.disc <- names(which(dum.disc == TRUE))

  # Check their concentartion
  disc.df <- data.frame(Concentration = sapply(tabular[dum.disc], function(x) {
    round(x[1] / sum(x), 3)
  }))
  
  # Find collinear variables
  coll.vars <- find.collinear(mm)
  
  ## Deal w/ originally dichotmous variables
  # Originally dichotmous variables will be discarded
  bin.disc <- grep("not mentioned", dum.disc, value = TRUE) # identify
  bin.disc <- gsub("not mentioned", "", bin.disc) # clean
  bin.disc <- c(bin.disc, "v169")
  EVS2017_final <- EVS2017_final[ , !colnames(EVS2017_final) %in% bin.disc] # get rid of
  
  ## Deal w/ categorical "non-variables"
  # Dummies for unpopular categories need to be dealt on a case by case basis
  
  # Religiosity
    round(prop.table(table(EVS2017_final$v51v52_comb)), 3)
    
    # Collapse the very unpopular categories
    denom <- fct_collapse(EVS2017_final$v51v52_comb,
                          protestant = c("Protestant"),
                          christian = c("Roman catholic"),
                          none = c("No Religion",
                                   "Orthodox",
                                   "Buddhist", 
                                   "Free church/Non-conformist/Evangelical",
                                   "Muslim",
                                   "Hindu",
                                   "Jew",
                                   "Other")
    )
    round(prop.table(table(denom)), 3)
    
    # Replace variable in original data
    EVS2017_final$v51v52_comb <- denom
    
  # v234 marital status
    prop.table(table(EVS2017_final$v234))
    v234 <- fct_collapse(EVS2017_final$v234,
                         partner = c("married",
                                     "registered partnership"),
                         never = c("never married and never registered partnership"),
                         past = c("divorced", "widowed", "separated")
    )
    prop.table(table(v234))
    EVS2017_final$v234 <- v234
    
  # v238: do you live with your parents/parents in law
    prop.table(table(EVS2017_final$v238))*100
    v238 <- fct_collapse(EVS2017_final$v238,
                         yes = c("yes, own parent(s)",
                                     "yes, both own parent(s) and parent(s) in law",
                                     "yes, parent(s) in law"),
                         no = c("no")
    )
    prop.table(table(v238))*100
    EVS2017_final$v238 <- v238
    
  # v244: 
    prop.table(table(EVS2017_final$v244))*100
    v244 <- fct_collapse(EVS2017_final$v244,
                         full = c("30h a week or more", "military service", "self employed"),
                         part = c("less then 30h a week"),
                         none = c("unemployed", "other", 
                                  "student", "disabled",
                                  "homemaker not otherwise employed")
    )
    prop.table(table(v244))*100
    EVS2017_final$v244 <- v244
    
  # v246: SES
    prop.table(table(EVS2017_final$v246))*100
    v246 <- fct_collapse(EVS2017_final$v246,
                         I = c("I :Higher Controllers"),
                         II = c("II :Lower Controllers"),
                         III = c("IIIa:Routine Nonmanual", 
                                 "IIIb:Lower Sales-Service"),
                         IV = c("IVa:Selfempl with empl",
                                "IVb:Selfempl no empl",
                                "IVc:Selfempl Farmer"),
                         V_VI = c("V :Manual Supervisors",
                                  "VI :Skilled Worker"),
                         VII = c("VIIa:Unskilled Worker", "VIIb:Farm Labor")
    )
    prop.table(table(v246))*100
    EVS2017_final$v246_egp <- v246
    
    dim(EVS2017_final)
    
# Save Results of cleaning and imputation ---------------------------------
    
  out <- list(orig = EVS2017_cl,
              full = EVS2017_final,
              note = paste0("2017 EVS dataset for selected countries and ",
                            "variables, ready for resampling study.",
                            " The two dataset have different dim because",
                            " after imputation a check for collinear vars",
                            " and constants is performed to exclude them") 
              )
  
  saveRDS(out, "../data/exp4_EVS2017_full.rds")
  
# Prepare Pattern For Missing Data ----------------------------------------
# According to Matrix Desing
  
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
  
  # Verify proportions
  table(var_class_f$Block)
  
  # Store results
  saveRDS(var_class_f, "../data/exp4_EVS2017_md.rds")

  