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
	* Variables supported: continuous, dichotomous, and categorical covaraites
		* IVs continuous, dichotomous and categorical (and mixed!);
		* DVs continuous, dichotomous and categorical;
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
	* Variables supported: 
		* IVs continuous and dichotomous (but not mixed!);
		* DVs continuous and dichotomous;
	* Limitations: 		
		* (1) requires the imputer to know the analysis model before imputation (look at the paper's 
			imputation model);
		* (2) requires a fully observed y variable. For the future, you could take insipiration 
			from the Sequential BART paper but avoid the assumption of knowing the analysis model beforehand;
		* (3) can only generate 5 multiply imputed datasets in current implementation.
	* Packages: A 'sbart' package was developed but has been removed from CRAN (and no updates since 2 years). 
		On github you have found 'bartpkg1' which is basically 'sbart'. The 'bart' function in BayesianTrees is 
		not directly applicabale as I cannot update the dataset that is used at each iteration of the mcmc
		algorithm. You could scrape the basic tree growing parts of bart from either backage and build 
		your own implementation of this algorithm.

* **MICE-Random Forest** following Shah et al 2014 (same as bart but the difference is how the 
	conditional mean is found)

* **Frequentist regularized-MICE** (directly and indirectly) following Deng et al 2016;
	There are no references to specific packages that implemented the methods but there are detailed 
	descriptions of how to make regularized multiple imputations.
	
	*Implementation status: DURR*
	* Imputation Model Variables supported: 
		* IVs: continuous, dichotomous, multinomial, ordinal, mixed;
		* DVs: continuous, dichotomous, multinomial, ordinal
	* Limitations:
		* DURR implementation at the moment gives you a terrible lasso selection for gaussian case. FIX on monday
		* number of iterations is not clearly specified in articles. Says to keep the last few datasets after convergence
		  but no clear measure of convergence is provided.
	* Packages: none

	*Implementation status: IURR*
	* Imputation Model Variables supported: 
		* IVs: continuous, dichotomous, multinomial, ordinal, mixed;
		* DVs: continuous, dichotomous, ordinal.
	* Limitations:
		* currently this apporach does not support multinomial outcome variables because the 
			paper does not specify how to do it and the way I perform a multinomial logistic regression 
			with lasso penality will provide different active sets (variable selection) for each baseline 
			category logit, meaning that I cannot use it for variable selection in a straightforward 
			way. In the logistic regression case, what you do is identify the predictors of the logit with
			the penalised method, then fit a logistic model using just those predictors and estimate parameters
			with a regular ML estimation. Finally, you use the parameters covariance matrix estiamted with such
			ML model to sample from a nomral distribution random values of the parameters to deifne a predictive
			model for the missing values. When using a regularized method for multinomial logistic, a log-linear
			model is used and each poisson count (how many observations in are in response category j of the dv, 
			given X) are modelled with different predictors (because of how the lasso penality works in the glmnet
			package). Hence, a unique active set is not identified and I cannot fit a regular multinomial logit,
			nor log linear model, to obtain a covariance matrix of parameter estiamtes that I would need to perform
			the missing data imputation.
		* number of iterations as for DURR
	* Packages: none

* **Bayesian LASSO** following Zhao Long 2016;
	* Imputation Model Variables supported: 
		* IVs: continuous, dichotomous, multinomial, ordinal, mixed;
		* DVs: continuous.
	* Limitations:
		* cannot impute dichotomous or polytmous variables at the moment
	* Packages: 'monomvn' R package implements the BLasso parameters selection part, I wrote code for implementing the
		MICE-like imputation algorithm.

* **PCA based auxiliary variables inclusion** following Howard et al 2015
	It is extendable to categorical and ordinal variables but I'm not familiar with PCA on these type
	of data so I need to look into it (Jolliffe 1986 discusses some ideas).
	
	*Implementation status:*
	* Imputation Model Variables supported: 
		* IVs: continuous, dichotomous, and categorical, mixed (see bottom notes for how);
		* DVs: continuous, dichotomous, and categorical DVs (you simply run mice on a
			low-dimensional dataset of DV and IVs + auxiliary)
	* Limitations: 
		* The initialization procedure advised for in the paper is plausible only for a low dimensinal
		  case ("A single Markov chain Monte Carlo (MCMC) imputation was used (including X and Y) as an 
		  intermediate step to acquire a complete data set on the auxiliary variables for the subsequent 
		  principal component analysis." p 289 under "Extracting Principal Components"). We might want to
		  use the form of imputation for PCA presented in the factominer R package. For now no intialization
		  is required for the auxiliary variables, so I prepare a dataset with PCs and then use it in MICE
		  directly. In the future, for each variable with missing values, a dataset with only the variable
		  and the first n PC components could be used (n should be defiend, look into the paper again) where all
		  other variables (interaction and poly terms) are used to compute the PC as well.
		* A priori knowledge of analysis model - Related to this issue, the method proposed by
		  Howard et al 2015 extracts PCs only form the auxiliary variables and not from the real or otherwise
		  selected predictors in the models and then uses the predictors + auxiliary Principal Components for
		  the actual imputation. This however, does not allow to automatically include interaction effects and
		  squared terms between real / selected predictors and auxiliary variables. In the implementation
		  I perfromed here I decided to extract PCs
		  form the auxiliary variables which have no missing values and then use the predictors + auxiliary PCs
		  for imputation with regular MICE. This is howver undesirable as we are not exlpoiting the full potential
		  of PCA for dimension reduction; this would still require to manually specify interactions and squared 
		  terms between real/selected predictors in the imputation model; we need to assume that we know the analysis 
		  model (ie which variables are auxiliary and which variables are predicotrs at least) before we 
		  perform the imputation, which is not ideal.
		* I have a bare bone implementation: no functions, just continuous, only squared terms but PcAux provides
		  a complete implementation that resolves the two aforementioned issues
	* Packages/Software: 
		* PcAux package implements Howards et al 2015
		* look into the missMDA function of the R package 'factominer'. It might be doing something close
		  to what you want (but it might also just be a form of MI *for*, rather than *using*, PCA).
	* Notes:
		* the single imputation referered to in PcAux package is one run of MICE with m = 1 and pmm method.
		* in PcAux the interaction variables missingness is done according to a transform then impute approach
		  (as Von Hippel 2009 advised for)
		* how are ordinal, non-ordinal, and dichotomous factors treated by PcAux? ordinal factors are casted to numeric
		  variables and treated as any oder continuous variable; nominal factors are dummy coded and just left as such
		  along with other potential dichotomous variables
		* in PcAux PCs are computed once, then used iteratively to predict the missing values in a MICE procedure
		  but their values stay constant across iterations and different data imputations.
		* number of PCs to be used: in the simualtion Howard et al 2015 use always just 1, in the real data 
		  example they use different numbers to see what is the effect of such choice and they conclude for example
		  that using 14 PCs that explain 55% of the auxiliary variables variance or 7 PCs that explain 40$ yields 
		  essentialy equivalent FMIs. In my implementation I will make the conservative choice of including a higher
		  number of PCs (say the ones that explain 50% of the variance, although why not include all?)
		* future perspective: I would like to have a version that uses all variables to compute the PC scores, 
		  this would help including in a automatic way the interaction and poly terms. This would require some
		  repetitions of PCA at least once per variable with missing values, if not once per variable per 
		  iteration (needs more thought); what about ICA?


