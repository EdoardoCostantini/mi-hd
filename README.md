# imputeHD-comp
Repository hosting project high-dimensional imputation comparison.

## Summary
Reasearchers working with large social surveys need tools to correct for the bias introduced by nonresponses.
The large number of items recorded in these surveys, coupled with their longitudinal nature and the necessity
of preserving complex interactions and nonlinear relations, easily produces high-dimensional ($p>n$) imputation 
problems that prevent a straightforward application of imputation algorithms such as MICE.

Furthermore, when employing Multiple Imputation to deal with missing values, data handlers tend to prefer including more
predictors in the imputation models as to reduce chances of uncongenial imputation and analysis models.
High-dimensional data imputation settings represent both an obstacle and an opportunity in this sense: 
an obstacle, as in the presence of high-dimensional data it is simply not possible to include all variables in 
standard parametric imputation models; 
an opportunity, because the large amount of features available has the potential to reduces the chances of 
leaving out of the imputation models important predictors of missignenss.

With this study we set out to present a thorough review and comparison of MI approaches for high-dimensional 
datasets. The goal is assessing how they meet the requirements of statistical validity of the analysis 
performed on the treated data.


