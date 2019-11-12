Here is a list of the methods we have decided to inlcude in the comparison. 
For each method I will soon write down the state of their implementation
for the specific goals of this project (e.g. what types of variables are 
supported, highdimensionality etc.)

* Sequential MI CART following Reiter 2005, and Drechsler Reiter 2011.

* BART and Bayesian CART imputation (if we want polytomous categorical) following Xu et al 2016.
	BART implementation w/ the bayesTree R package is functioning for both coninuous and binary 
	target variables and highdimensional settings.

* MICE-Random Forest following Shah et al 2014 (same as bart but the difference is how the 
	conditional mean is found)

* Frequentist regularized-MICE (directly and indirectly) following Chang et al 2016;

* Bayesian LASSO following Zhao Long 2016;

* PCA based auxiliary variables inclusion following Howard et al 2015


