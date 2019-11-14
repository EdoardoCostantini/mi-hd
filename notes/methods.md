Here is a list of the methods we have decided to inlcude in the comparison. 
For each method I will soon write down the state of their implementation
for the specific goals of this project (e.g. what types of variables are 
supported, highdimensionality etc.)

* Sequential MI CART following Reiter 2005, and Drechsler Reiter 2011.
	I can use any CART growing R package. Supports any type of outcome variable without any effort

* BART and Bayesian CART imputation (if we want polytomous categorical) following Xu et al 2016.
	BART implementation w/ the bayesTree R package is functional for both coninuous and binary 
	target variables and highdimensional settings.
	With BART, you obtain an empirical conditional distribution of a Y given some set of 
	covariates X. If we grow 1 tree T, we obtain b terminal nodes, and in each a value mu defines
	the conditional mean of Y given a precise value of X. If we grow a set of m trees (say 10) we
	obtain 10 different expected values of Y for each value of X. The sum of such different 
	conditional means is the best estimate we have of the conditional mean of y given a specific value
	of X. These m trees are sampled from a posterior distribution that defines the porbability of each 
	of m fixed number of trees given the data we have observed. With priors we can keep the individual 
	effects of each of m trees small. The MCMC algorithm that allows to sample these m trees makes you
	sample each tree conditionally on all the other trees. After a burnin period, we may consider a set of
	m trees obtained with 1 iteration as a random sample form the posterior distribution of the sets of m
	trees fiven the data observed. Each of the K sets of m trees obtained in the after burn-in iterations
	provides one estimate of the conditional mean of y given X. We can then approximate the E(f(x)|y). 
	Each iteration provides a set of m trees, and a draw from the empirical distribution of f(x)|y.

	Using the 'bart' function with y.train = variable with missing values, x.train = all other 
	variables, we obtain ndpost/keepevery draws from the posterior distribution of the missing 
	values given observed part of y.train, all other vairables, and model specifications 
	(in yhat.train, and also in yhat.test if you use the same X for x.train and x.test)

* MICE-Random Forest following Shah et al 2014 (same as bart but the difference is how the 
	conditional mean is found)

* Frequentist regularized-MICE (directly and indirectly) following Deng et al 2016;
	There are no refere3nces to specific packages that implemented the methods but there are detailed 
	descriptions of how to make regularized multiple imputations.

* Bayesian LASSO following Zhao Long 2016;

* PCA based auxiliary variables inclusion following Howard et al 2015
	Supports auxiliary vairables with missing values
	It is extendable to categorical and ordinal variables but I'm not familiar with PCA on these type
	of data so I need to look into it (Jolliffe 1986 discusses some ideas).


