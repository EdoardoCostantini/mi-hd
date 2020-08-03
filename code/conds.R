### Title:    Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-05-19
### Notes:    some ideas on the conditions that I could specify for
###           the simulation study
# Simulation characteristics 
# Fixed and vairable factors

# Data generation ---------------------------------------------------------

# Sample size
n <- c(200)
  # fixed seems a no brainer to me, but how big?
  # should mimic an EVS study? I think not, as the
  # EVS resampling study will do that for us.
# Number of predictors
p <- c(50, 200, 500)
  # well behaved n > p, up to n << p
  # this is the main factor. How do the methods behave as the
  # the proportion of important predictors decreases.
  # e.g. expect PCA to perform worse and worse, Blasso and other reg
  # to stay put, random forest to be impossible, 
  # This way we can address better the point of psychology data1

# Generation of X
atuoRho <- c(.5)
# it is good to think about collinearity bettter to think of a latent structure
# generate X accirdubg ti sine kind of laatent structure (in gorups) and manipulate 
# how important the groupings are and generate Y based on the latent vairables and not on the 
# X (maybe look at Elasti net paper for orignal idea but the the latent vairables there
# are basically perfect measurements).
  # in general, this will impact the correlation between 
  # auxiliary and imputation outcome variables.
  # Why would I vary this factor as in Zhao Long 2016?
X_vars <- runif(p, 0, 10)
X_covs <- runif(p*p, 0, 1)

# True substantive data generation
# Substantive model interaction fixed to yes
y1 = b0 + b %*% X[, 1:7] + b_i1 * X[,1]*X[,2]
y1 = b0 + b %*% X[, 1:7] + b_p1 * X[,1]**2
  # Interesting to have an interaction coefficient in the analysis model
  # to see how imputation methods influences performance for this estimate.

# N var w/ miss
# How many variables will have missing values. We want to focus on
# multivariate missingness. More then one variables w/ missingness
# is a fixed aspect, but which and how many needs thinking.
# 1. Analysis model(s) DV
# 2. Analysis model(s) IVs
# 3. Auxiliary variables
# At the moment, I would fix this to be: all of variables in 1 and 2, 
# plus some in 3
pm <- c(.1, .2, .3) # proportion of missing cases because in sparses situation some could perfoem better
N_varmis <- c(.5)
# Y variable included as missing variable
# AS_miss_vars  (q in ZhaoLong2016)
# How many variables are used in the generation of the 
# variables that will have missing values
# Category B variables from Collins Et Al 2001 paper.
#  correlated with Y and is not a correlate of missingness
# Interest: 
# 1. If included in imputation models, reduce uncertainty
#    in the data -> smaller standard errors, greater power,
#    narrower CI w/out hurting coverage
#    [see Collins Et Al 2001, p. 338]
# 2. Even MI-TRUE does not perform well if the data gen
#    includes many category b variables (large bias). 
#    Blasso, DURR and IURR mantain low bias and 
#    good coverage! (Zhao Long 2016, table 1 and 3)

# 1. Number of variables influnecing scores
AS_varmis <- c(5, 50)

# 2. Complexity of relationship
# INteraction terms
# No/yes interactions in the generating model for the variables
# that will have missing values
# This is fundamental especially beacuse the tree models have this as
# the perk. And they have to be nontrivial
# And you would have an imputation set up with raw data and one 
# with no 
INT_varmis <- c(TRUE, FALSE)
# Polynomial terms
POL_varmis <- c(TRUE, FALSE) 

# THink about sparsity of the data generating model
# And fix P and iteractions when you are looking at this
# and cross q sizes for data generation and missing data imposition
# Maybe just two levels sparse and not sparse in Baysesian EN paper.
# Good to see the missing data genartion Active Set.

# Think about it more as in experiments around a theme.
# You can do it in a similar way because simulations are just 
# factorial experiments,
# Figure out the least interesting crossings and marginalized then out.

# Missing data mechanism --------------------------------------------------

# The advantage of including many auxiliary variables
# is expected to be important the more complex the missingness 
# generating model is. This complexity depends on:
# 1. Number of variables influencing missingness
  AS_prbmis <- c(3, 5, 10)
# 2. Complexity of relationships (i.e. interaction and polytomous 
#    terms). This could be expressed in different ways:
# Include as a cross condtion to know if there is a benefit in considering 
# it or not
# How methods deals with undetermined sismtes. Does not need to be crossed with the interaction question 
# The question of recovering interactions is 
# With a fixed p you corss the codntions of where the interaction is.
  # A) degree of interactions / poly terms
  INT_prbmis <- c(0, .5, 1)    # proportion of variables interacting
  POL_prbmis <- c(0, 1, 2)     # number of variables with poly terms
  # B) presence / absence of interactions and polu terms
  INT_prbmis <- c(TRUE, FALSE) 
  POL_prbmis <- c(TRUE, FALSE)
  # C) general complexity of missingness (following Collins et al 2001)
  MAR <- c("LINEAR", "INTER") # type of true relatioship
  
# Whether or not data should be augmented w/ interactions and polyterms
# before imputation
  aug <- c("TRUE/FALSE")

# Planned missingens
  plannedMCAR <- c(FALSE, TRUE)
  
# Latent Structure
  latent <- c(FALSE, TRUE)
  
# Experiment 1.1 ----------------------------------------------------------
# Name:     Increase in dimensionality
# Question: Which data management / method combination allows for
#           best recovery of interaction terms?
# Notes     1. Do the methods perform well in low dimensional cases?
#             (compare with traditional mice)
#           2. Which methods keep stronger performances while dimesionality
#             and proportion of missings increase?
  
  expand.grid(n = n, 
              latent = latent[1],
              atuoRho = atuoRho,
              MAR = MAR[1],
              MCAR = plannedMCAR[1],
              p = p,
              pm = pm)
  
# Experiment 1.2 ----------------------------------------------------------
# Name:     Increase in dimensionality
# Question: Which data management / method combination allows for
#           best recovery of interaction terms?
# Notes     1. Do the methods perform well in low dimensional cases?
#             (compare with traditional mice)
#           2. Which methods keep stronger performances while dimesionality
#             and proportion of missings increase?
  
  expand.grid(n = n, 
              latent = latent[1],
              MAR = MAR,
              p = p,
              pm = .2)
  
# Experiment 2 ------------------------------------------------------------
# Name:     Recovering Interaction relationshps: 
# Question: Which data management / method combination allows for
#           best recovery of interaction terms?
# Notes:    1. You want to see how the methods deal with no data augmentation
#             with interactions present in different points of the data gen
#             Ideally, things like CART and RF will pick up that relationship
#             while the other will not (and this will impact their performances)
#           2. You want to see that Random Forest and CART cannot proceede well
#             in data augmentation scenarios.
#           3. You want to compare the methods at their best, e.g. IURR with
#             data augmentation, and Random Forest without DA.

  expand.grid(n = n, 
              p = p[2],
              pm = pm[3], 
              atuoRho = atuoRho,
              MAR = MAR[1],
              AS_varmis = AS_varmis[1], 
              AS_prbmis = AS_prbmis[2], 
              INT_varmis = INT_varmis,
              POL_varmis = POL_varmis[2],
              INT_prbmis = INT_prbmis,
              POL_prbmis = POL_prbmis[2],
              augment = aug)

# Experiment 2.b ------------------------------------------------------------
# Name:     Recovering Polynomial relationshps: 
# Question: Which data management / method pair 
#           allows for best recovery of polynomial terms?
# Note:     1. You  want to see basically the same 3 points as in Exp 2
#             but with polynomials instead of interactions.
#           2. This allows you to check non-linear MAR, so you do not have to
#             to add the linear/convex mar in Experiment 1. Of course here
#             you are not seeing how it interplays with data dimensinality
#             or proportion of missingness, but it's enough to see how 
#             linear / non-linear mar is managed in a difficult set up
#             (say n = p, pm = .3)
  
  expand.grid(n = n, 
              p = p[2],
              pm = pm[3], 
              atuoRho = atuoRho,
              MAR = MAR[1],
              AS_varmis = AS_varmis[1], 
              AS_prbmis = AS_prbmis[2], 
              INT_varmis = INT_varmis[2],
              POL_varmis = POL_varmis,
              INT_prbmis = INT_prbmis[2],
              POL_prbmis = POL_prbmis,
              augment = aug)

    
## Option 2 -  More generic interaction factor
intType <- c("NONE", "SUB", "VARMIS", "PRBMIS", "ALL")
polyType <- c("SUB", "VARMIS", "PRBMIS", "ALL") # second order

expand.grid(n = n, 
            N_varmis = N_varmis, 
            MECH = MECH,
            atuoRho = atuoRho,
            AS_varmis = AS_varmis, 
            AS_prbmis = AS_prbmis, 
            intType = intType,
            p = p)

## Option 3 - Linear non linear mar is more important

MAR <- c("LINEAR", "INTER", "CONVEX")

expand.grid(n = n, 
            N_varmis = N_varmis, 
            atuoRho = atuoRho,
            AS_varmis = AS_varmis, 
            AS_prbmis = AS_prbmis, 
            intType = intType[c(-4)],
            p = p,
            MECH = MAR)
