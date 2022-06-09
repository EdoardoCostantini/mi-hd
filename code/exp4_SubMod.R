### Title:    Compare dataset with missings and pmm imputed dataset
### Author:   Anonymized for peer review
### Created:  2020-09-18

rm(list=ls())
source("./init_general.R")
source("./exp4_init.R")

library(ggplot2)

# Data Prep ---------------------------------------------------------------

EVS_dt <- readRDS("../data/exp4_EVS2017_full.rds")

dt.o <- EVS_dt$orig
dt.f <- EVS_dt$full

# Distributional Checks ---------------------------------------------------
  par(mfrow=c(1,2))

  # Age
  lapply(list(dt.o, dt.f), 
         function(x) plot(density(x$age,
                                  adjust = 2,
                                  na.rm = TRUE),
                          xlim = c(0, 100), xlab = "Age")
  )
  
  # Left Right party: v174_LR
  lapply(list(dt.o, dt.f), function(x) plot(density(x$v174_LR,
                                                    adjust = 2,
                                                    na.rm = TRUE),
                                            xlim = c(0, 11))
  )
  
  # SES: v246_egp
  lapply(list(dt.o, dt.f), function(x) plot(density(x$v246_egp, 
                                                    adjust = 3,
                                                    na.rm = TRUE),
                                            xlim = c(0, 20)
  )
  )
  
  # INCOME: v261_ppp
  lapply(list(dt.o, dt.f), function(x) plot(density(x$v261_ppp, 
                                                    adjust = 5,
                                                    na.rm = TRUE),
                                            xlim = c(0, 7)))
  
  summary(lm("v174_LR ~ v261_ppp + v246_ISEI +", data = dt.f))
  
  summary(lm("v174_LR ~ v261_ppp + v246_ISEI + age", data = dt.f))
  summary(lm("v6 ~ v153 + v154 + v155", data = dt.f))
  
  summary(
    lm("v6 ~ v153 + v154 + v155", 
       data = as.data.frame( scale(dt.f[, c("v6", "v153", "v154", "v155")]) )
    )
  )

# Euthanasia --------------------------------------------------------------
# Inspired by Koneke 2014 (model 2)
  # Euthanasia
  table(dt.f$v156)
  lapply(list(dt.o, dt.f), 
         function(x) plot(density(x$v156,
                                  adjust = 5,
                                  na.rm = TRUE),
                          xlim = range(dt.f$v156)
         )
  )
  euth <- dt.f$v156
  
  # Country
  country <- dt.f$country
    
  # General Trust
  table(dt.f$v31)
  trust.g <- dt.f$v31
  
  # Confidence in Health care sys
  table(dt.f$v126) 
  lapply(list(dt.o, dt.f), 
         function(x) plot(density(x$v126,
                                  adjust = 5,
                                  na.rm = TRUE),
                          xlim = range(dt.f$v126)
                          )
         )
  
  trust.hs <- match(dt.f$v126, 4:1)
  cbind(dt.f$v126, trust.hs)
  trust.hs.o <- dt.f$v126
  
  # Confidence in press
  table(dt.f$v118)
  lapply(list(dt.o, dt.f), 
         function(x) plot(density(x$v118,
                                  adjust = 5,
                                  na.rm = TRUE),
                          xlim = range(dt.f$v118)
         )
  )
  trust.pr <- match(dt.f$v118, 4:1)
  cbind(dt.f$v118, trust.pr)
  trust.pr.o <- dt.f$v118
  
  # Confidence in state scale (need to create scale)
  # parlianment, police, justice sys, goverment
  # Distributions are ok
  dat.trust <- lapply(list(dt.o=dt.o, dt.f=dt.f), 
         function(x) {
           x$trust.state <- rowMeans(x[, c("v121", "v120", "v127", "v131")])
           return(x)
         }
  )
  
  lapply(dat.trust, 
         function(x) plot(density(x$trust.state,
                                  adjust = 5,
                                  na.rm = TRUE),
                          xlim = range(dat.trust$dt.f$trust.state)
         )
  )
  
  # But recode in more meaningful way
  trust.s.var <- c("v121", "v120", "v127", "v131")
  trust.s.o <- rowMeans(dt.f[, trust.s.var]) # original scaling
  invCod <- sapply(trust.s.var, 
                   function(x) match(dt.f[, x], 4:1))
  dt.f[, trust.s.var] <- invCod
  
  trust.s <- rowMeans(dt.f[, trust.s.var]) # inverted scaling
  
  # Covariates
  # Age
  lapply(list(dt.o, dt.f), 
         function(x) plot(density(x$age,
                                  adjust = 2,
                                  na.rm = TRUE),
                          xlim = c(0, 100), xlab = "Age")
  ) 
  
  age <- dt.f$age
  
  # Education
  table(dt.f$v243_ISCED_1)
  lapply(list(dt.o, dt.f), 
         function(x) plot(density(x$v243_ISCED_1, 
                                  adjust = 5,
                                  na.rm = TRUE),
                          xlim = range(dt.f$v243_ISCED_1))
  )
  edu <- dt.f$v243_ISCED_1
  
  # Gender (1 = male)
  sex <- dt.f$v225
  
  # Religiousness
  table(dt.f$v6)
  cbind(dt.f$v6,
        match(dt.f$v6, 4:1)-1)
  rel <- match(dt.f$v6, 4:1)-1
  rel.o <- dt.f$v6
  
  # Religious Denomination
  denom <- dt.f$v51v52_comb
  
# > Models ####
  
  # Model 1
  mod1 <- lm(euth ~  trust.g + trust.hs + trust.pr + trust.s + country)

  round(summary(mod1)$adj.r.squared*100, 3)
  round(summary(mod1)$coefficients, 3)[, 1:2]
  
  # Model 2
  mod2 <- lm(euth ~
               trust.g + trust.hs + trust.pr + trust.s + edu + 
               age + sex + rel + denom + country)
  round(summary(mod2)$adj.r.squared*100, 1)
  round(summary(mod2)$coefficients, 3)[, c(1:2, 4)]
  
  mod2 <- exp4_fit_mod1(list(GS = dt.f))$GS
  mod2
  round(summary(mod2)$adj.r.squared*100, 1)
  round(summary(mod2)$coefficients, 3)[, c(1:2, 4)]

# Left Right Vote Gender --------------------------------------------------
# Inspired by ImmerzeelEtAl2016
  
# Left / Right voting
  table(dt.f$v174_LR)
  lapply(list(dt.o, dt.f), 
         function(x) plot(density(x$v174_LR, 
                                  adjust = 5,
                                  na.rm = TRUE),
                          xlim = range(dt.f$v174_LR))
  )
  lr <- dt.f$v174_LR
  
# Country
  country <- dt.f$country  
  
# Female
  female <- relevel(dt.f$v225, ref = "male")
  
# Employment Status  
  # employed, self-employed and inactive/unemployed.
# Occupation
  SES <- dt.f$v246_egp

# Native attitudes (mean of itmes)  
  nativ <- c("v185", # jobs
             "v186", # crime
             "v187" # strain on welfare
  )
  NatAt <- rowMeans(dt.f[, nativ])
# Authoritarian Attitudes

# Low and order attitudes
  strongL <- match(dt.f$v145, 4:1)
  order <- forcats::fct_collapse(dt.f$v110,
                               no = c("fighting rising prices", 
                                       "more say in important government decisions", 
                                       "protect freedom of speech"),
                               yes = c("maintaining order in nation")
                               )
  
# Political Interest
  polInterest <- match(dt.f$v97, 4:1)
  
# Political Action
  action <- paste0("v", 98:101)
  polAction <- rowMeans(dt.f[, action])
  polAction_r <- case_when(
    polAction %in% 1 ~ 3, # have done all
    polAction %in% seq(1.25, 2.75, by = .25) ~ 2,
    polAction %in% 3 ~ 1 # would never do any
  )
  cbind(polAction, polAction_r)
  
# Covariates
  age <- dt.f$age
  edu <- dt.f$v243_ISCED_1
  mat <- forcats::fct_collapse(dt.f$v234,
                               partner = c("married", "registered partnership"),
                               never = "never married and never registered partnership",
                               had = c("divorced", "separated", "widowed"))
  
  urb <- dt.f$v276_r # dummy coded
  
  rel <- factor(case_when(
    dt.f$v54 %in% 7 ~ 0, # less than once a year
    dt.f$v54 %in% 4:6 ~ 1, # sometimes
    dt.f$v54 %in% 1:3 ~ 2, # once per month or more
  ), labels = c("never", "sometimes", "often"))
  
  denom <- dt.f$v51v52_comb
  
# > Models ####
  
  mod3 <- lm(lr ~ female + SES + NatAt + strongL + order + polInterest + 
       polAction_r + age + edu + mat + urb + denom + country)
  
  summary(mod3)
  round(summary(mod3)$adj.r.squared*100, 1)
  round(summary(mod3)$coefficients, 3)[, c(1:2, 4)]

  mod3 <- exp4_fit_mod2(list(GS = dt.f))$GS
  coef(mod3)
  confint(mod3)
  round(summary(mod3)$adj.r.squared*100, 1)
  round(summary(mod3)$coefficients, 3)[, c(1:2, 4)]
  
  