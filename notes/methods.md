**Methods to be compared**

Here is a list of the methods we have decided to inlcude in the comparison. 
For each method I will soon write down the state of their implementation
for the specific goals of this project (e.g. what types of variables are 
supported, highdimensionality etc.)

* **Sequential MI CART** following Reiter 2005, and Drechsler Reiter 2011.
	I can use any CART growing R package. Supports any type of outcome variable without any effort
	There is no guideline on what initialization method to use as most of the papers that proposed this method 
	deal w/ creating synthetic datasets, where all variables are actually observed.
	
	*Implementation status:*
	* Variables supported: continuous, dichtomous, and categorical covaraites
	* Limitations: none
	* Software/Packages: 'mice' basic mice imputation function with a 'cart' option does everything you need; 
	you have also implemented a version of it that uses bayesian bootstrap.

* **BART** and Bayesian CART imputation (if we want polytomous categorical) following Xu et al 2016.
	BART implementation w/ the bayesTree R package is functional for both coninuous and binary 
	target variables and highdimensional settings.
	With BART, you obtain an empirical conditional distribution of a Y given some set of 
	covariates X. If we grow 1 tree T, we obtain *b* terminal nodes, and in each a value $\mu$ defines
	the conditional mean of Y given a precise value of X. If we grow a set of *m* trees (say 10) we
	obtain 10 different expected values of Y for each value of X. The sum of such different 
	conditional means is the best estimate we have of the conditional mean of y given a specific value
	of X. These *m* trees are sampled from a posterior distribution that defines the porbability of each 
	of *m* fixed number of trees given the data we have observed. With priors we can keep the individual 
	effects of each of m trees small. The MCMC algorithm that allows to sample these *m* trees makes you
	sample each tree conditionally on all the other trees. After a burnin period, we may consider a set of
	m trees obtained with 1 iteration as a random sample form the posterior distribution of the sets of m
	trees fiven the data observed. To predict the Y value for one particular $x_i$, we can take the mean 
	of the sums of $\mu$ obtained in each of the K sets of *m* trees sampled.

 	*Implementation status:*
	* Variables supported: continuous, dichtomous, and categorical covaraites
	* Limitations: (1) requires the imputer to know the analysis model before imputation. (2) requires 
	a fully observed y variable. For the future, you could take insipiration from the Sequential BART paper 
	but avoid the assumption of knowing the analysis model beforehand.
	* Packages: A 'sbart' package was developed but has been removed from CRAN (and no updates since 2 years).
	The 'bart' function in BayesianTrees is not directly applicabale as I cannot update the dataset that
	is used at each iteration of the mcmc algorithm. You could scrape the basic tree growing parts of bart 
	from either backage and build your own implementation of this algorithm.

* **MICE-Random Forest** following Shah et al 2014 (same as bart but the difference is how the 
	conditional mean is found)

* **Frequentist regularized-MICE** (directly and indirectly) following Deng et al 2016;
	There are no refere3nces to specific packages that implemented the methods but there are detailed 
	descriptions of how to make regularized multiple imputations.

* **Bayesian LASSO** following Zhao Long 2016;

* **PCA based auxiliary variables inclusion** following Howard et al 2015
	Supports auxiliary vairables with missing values
	It is extendable to categorical and ordinal variables but I'm not familiar with PCA on these type
	of data so I need to look into it (Jolliffe 1986 discusses some ideas).
	
	*Implementation status:*
	* Variables supported: continuous, dichtomous (?), and categorical covaraites (?)
	* Limitations: NA
	* Packages: look into the missMDA function of the R package 'factoineR'. It might be doing something close
		to what you want (but it might also just be a form of MI *for*, rather than *using*, PCA)


