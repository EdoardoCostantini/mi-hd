### Title:    Replication Howard et al 2015
### Author:   Edoardo Costantini
### Created:  2020-05-18

# Lavan models
# NO AUXILIARY VARIABLES
# Create descriptive model object
mod_NAV <- '
    # means
    Y ~ 1
    X ~ 1
    
    # variances
    Y ~~ Y
    X ~~ X
    
    # covariances/correlations
    Y ~~ X
    '
# INCLUSIVE STRATEGY
mod_IAV <- '
    # means
    Y ~ 1
    X ~ 1
    A1 ~ 1
    A2 ~ 1
    A3 ~ 1
    A4 ~ 1
    A5 ~ 1
    A6 ~ 1
    A7 ~ 1
    A8 ~ 1
    
    # variances
    Y ~~ Y
    X ~~ X
    A1 ~~ A1 
    A2 ~~ A2 
    A3 ~~ A3 
    A4 ~~ A4 
    A5 ~~ A5 
    A6 ~~ A6 
    A7 ~~ A7 
    A8 ~~ A8 
    
    # covariances/correlations
    Y ~~ X + A1 + A2 + A3 + A4 + A5 + A6 + A7 + A8
    X ~~ A1 + A2 + A3 + A4 + A5 + A6 + A7 + A8
    A1 ~~ A2 + A3 + A4 + A5 + A6 + A7 + A8
    A2 ~~ A3 + A4 + A5 + A6 + A7 + A8
    A3 ~~ A4 + A5 + A6 + A7 + A8
    A4 ~~ A5 + A6 + A7 + A8
    A5 ~~ A6 + A7 + A8
    A6 ~~ A7 + A8
    A7 ~~ A8
    '
# PCA STRATEGY
mod_PCA <- '
    # means
    Y ~ 1
    X ~ 1
    PC1 ~ 1
    
    # variances
    Y ~~ Y
    X ~~ X
    PC1 ~~ PC1
    
    # covariances/correlations
    Y ~~ X + PC1
    X ~~ PC1
    '
