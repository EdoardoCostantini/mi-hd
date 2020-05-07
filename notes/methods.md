**Methods to be compared**

Here is a list of the methods we have decided to inlcude in the comparison. 
For each method I will soon write down the state of their implementation
for the specific goals of this project (e.g. what types of variables are 
supported, highdimensionality etc.)

* **Sequential MI CART (DONE)** following Reiter 2005, and Drechsler Reiter 2011.
	I can use any CART growing R package. Supports any type of outcome variable without any effort
	There is no guideline on what initialization method to use as most of the papers that proposed this method 
	deal w/ creating synthetic datasets, where all variables are actually observed.
	
	*Implementation status:*
	* Variables supported: continuous, dichotomous, and categorical covaraites
		* IVs continuous, dichotomous, categorical, ordinal as numeric (and mixed!);
		* DVs continuous, dichotomous, categorical, ordinal as numeric;
	* Limitations: none
	* Software/Packages: 'mice' basic mice imputation function with a 'cart' option does everything you need; 
	you have also implemented a version of it that uses bayesian bootstrap. You want to include your version 
	in the mice package. You have created  the mice.impute.cart.bb.R script in the miceImpHDv package version
	that is on your computer in the R-packages folder (under projects). Now you need to do the same you did 
	for the Random forests.

* **Frequentist regularized-MICE (DONE)** (directly and indirectly) following Deng et al 2016;
	There are no references to specific packages that implemented the methods but there are detailed 
	descriptions of how to make regularized multiple imputations.
	
	*Implementation status: DURR*
	* Imputation Model Variables supported: 
		* IVs: continuous, dichotomous, multinomial, ordinal, mixed;
		* DVs: continuous, dichotomous, multinomial, ordinal
	* Limitations:
		* number of iterations is not clearly specified in articles. Says to keep the last few datasets after convergence
		  but no clear measure of convergence is provided.
		* one good aspect of EN is that if Lasso or ridge perform better, then they are preferred; howver, the corssvalidation
		  of two parameters required by EN is extremely time consuming
	* Packages: none; R function to perform imputation is operational
	* Notes: Implemented the Elastic net exclusively because it reduces to LASSO or RIDGE that performs better. Theoretically,
		 EN is better at selecting variables from a group of variables that are highly correlated. As we are thinking of 
		 surveys as applications, EN seems theoretically more apt to the task beacuse of its grouping effect: if there are 
		 predictors that are highly correlated, then EN would tend to push their coefficient to be the same, while LASSO 
		 would tend to keep just one of them.
	* Doubts: Iterations and saving datasets for multivaraite missingness. I currently implemented this with parallel chains. 
		However, Deng et al 2016 seem to suggest that there is only 1 chain that converges after so many iterations and 
		that the last M dataset can be kept as multiple imputed datasets. Isn't there dependency that we want to avoid in 
		this method?

	*Implementation status: IURR*
	* Imputation Model Variables supported: 
		* IVs: continuous, dichotomous, multinomial, ordinal, mixed;
		* DVs: continuous, dichotomous, multinomial, ordinal.
	* Limitations:
		* number of iterations as for DURR
		* find an efficient crossvalidation for elastic net
		* weakness: when using this method it might happen that p > n variables are selected
			making the apporach fail.
	* Packages: none, R function to perform imputation is operational
	* Doubts: The sampling of sigma2. in your code there is no variation in the parameter sigma_j^2. You keep the 
		value estiamted with the lasso model. Ask Kyle for help.


* **Bayesian LASSO** following Zhao Long 2016;
	* Imputation Model Variables supported: 
		* IVs: continuous, dichotomous, multinomial, ordinal, mixed;
		* DVs: continuous.
	* Limitations:
		* cannot impute dichotomous or polytomous variables at the moment: needs a different algorithm 
	* Packages: 'monomvn' R package implements the BLasso parameters selection part, I wrote code for implementing the
		MICE-like imputation algorithm.
	* Notes: 
		* The lasso penalty is proportional to the log density of a laplacian distribution (or double exponential) (see 
		Tibshirani1996 section 5) and this distribution has a shape with mass in the center and in the tails meaning
		that when used as a prior for regression coefficeints, it will favor coefficients that are either 0 or large. 
		Advantages of this procedure is that it provides directly a valid posterior distribution of the coefficients 
		AND the posterior predictive distribution of the missing values.
		* HTLR R package has a very nice implementation of the bayesian multinomial logistic regression with 
		hyper-lasso prior for feature selection. Include it in the BLasso code for imputation and then you
		are good to go w/ it.

* **PCA based auxiliary variables inclusion (DONE)** following Howard et al 2015
	
	*Implementation status:*
	* Imputation Model Variables supported: 
		* IVs: continuous, dichotomous, and categorical, mixed (see bottom notes for how);
		* DVs: continuous, dichotomous, and categorical DVs (you simply run mice on a
			low-dimensional dataset of DV and IVs + auxiliary)
	* Limitations:  
		* A priori knowledge of analysis model - Related to this issue, the method proposed by
		  Howard et al 2015 extracts PCs only form the auxiliary vairables and not from the real or otherwise
		  selected predictors in the models and then uses the predictors + auxiliary Principal Components for
		  the actual imputation. This however, does not allow to automatically include interaction effects and
		  squared terms between real / selected predictors and auxiliary vairables. In the implementation
		  I perfromed here I decided to extract PCs
		  form the auxiliary variables which have no missing values and then use the predictors + auxiliary PCs
		  for imputation with regular MICE. This is howver undesirable as we are not exlpoiting the full potential
		  of PCA for dimension reduction; this would still require to manually specify interactions and squared 
		  terms between real/selected predictors in the imputation model; we need to assume that we know the analysis 
		  model (ie which vairables are auxiliary and which variables are predicotrs at least) before we 
		  perform the imputation, which is not ideal.
	* Packages/Software: 
		* PcAux package implements Howards et al 2015
		* look into the missMDA function of the R package 'factominer'. It might be doing something close
		  to what you want (but it might also just be a form of MI *for*, rather than *using*, PCA).
	* Notes:
		* the single imputation referered to in PcAux package is one run of MICE with m = 1 and pmm method.
		* in PcAux the interaction variables missingness is done according to a transform then impute approach
		  (as Von Hippel 2009 advised for)
		* how are ordinal, non-ordinal, and dichotomous factors treated by PcAux? ordinal factors are casted to numeric
		  variables and treated as any oder continuous variable if not specified otherwise with the ordVar argument?; 
		  nominal factors are dummy coded and just left as such along with other potential dichotomous variables (see 
		  dummyVars and probNoms arguments in the createPcAux function)
		* in PcAux PCs are computed once, then used iteratively to predict the missing values in a MICE procedure
		  but their values stay constant across iterations and different data imputations.
		* number of PCs to be used: in the simualtion Howard et al 2015 use always just 1, in the real data 
		  example they use different numbers to see what is the effect of such choice and they conclude for example
		  that using 14 PCs that explain 55% of the auxiliary vairable variance or 7 PCs that explain 40$ yields 
		  essentialy equivalent FMIs. In my implementation I will make the conservative choice of including a higher
		  number of PCs (say the ones that explain 50% of the vairance, although why not include all?)
		* future perspective: I would like to have a version that uses all variables to compute the PC scores, 
		  this would help including in a automatic way the interaction and poly terms. This would require some
		  repetitions of PCA at least once per variable with missing values, if not once per variable per 
		  iteration (needs more thought); what about ICA?

* **MICE-Random Forest (DONE)** following Shah et al 2014 (same as bart but the difference is how the 
	conditional mean is found) and DooveEtAl 2014

	*Implementation status:*
	* Variables supported: 
		* IVs everything;
		* DVs continuous and nominal (ordinal as continuous);
	* Limitations:
		* Ordinal variables are force to numericcaly continuous ones
	* Packages: (1) mice::mice.impute.rf() based on Doove et al 2014; (2) a modified version of mice.impute.rf() 
			from the mice package based on Shah et al 2014, available in my forked version of the mice package.
	* Notes: uses random forests to define the mean of the predictive distribution of the missing values and uses the
		out of bag mean squared error as variance. It differs from the regular mice::mice.impute.rf() in the following 
		way: The main difference is in how the multiple trees in the forest are used to impute the missing values: 
		in Shah et al’s algorithm for continuous variables, the missing value is sampled from a normal distribution 
		centred around the aggregated predicted value for the ’test’ data (imputation model X values for the unobserved 
		Y values) and the out-of-bag error is used as variance, while their categorical variable algorithm samples one 
		tree at random from the k trees that populate the forest and use it for prediction; in the ’mice’ package version, 
		each value of y_miss_j is sampled from a donor pool deﬁned by all of the y_obs_j that end up in the k leafs 
		deﬁned by the k trees in the random forest.
		The advantage comes from:
		* ensamble learning algorithm (not relying on a single tree, which is a greedy approach)
		* prediction improvement: reduces prediction variance, but also prediction bias

* **MissForest (DONE)**
*Implementation status:*
	* Variables supported: 
		* IVs everything;
		* DVs continuous and nominal (ordinal as continuous? not sure);
	* Limitations:
		* Single imputation
	* Packages: missForest R package 


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
		* (4) **UNDERSTANDING** - I still do not understand it fully. The main problem is what goes in the acceptance
			probability in step 1.b of appendix 2. If you solve that then maybe you can understand. You are close
			to understainding the code but still something is not there. Open the project in CLion and go to the
			main programme file: "mi-bart-y1", around line 200 the definition of the proposition acceptance
			probability numerator and denominator is done.
	* Packages: A 'sbart' package was developed but has been removed from CRAN (and no updates since 2 years). 
		On github you have found 'bartpkg1' which is basically 'sbart'. The 'bart' function in BayesianTrees is 
		not directly applicabale as I cannot update the dataset that is used at each iteration of the mcmc
		algorithm. You could scrape the basic tree growing parts of bart from either backage and build 
		your own implementation of this algorithm.
	* Notes: 

* **BART** my way.
	Given a basic model, Y = f(x) + E, with E ~ N(0, sigma_E^2), we can model f(x) as the sum of many trees by fitting
	a bart model (for Heteroskedastic version of BART see the rbart vignette). Fitting a bart model is done by running 
	a MCMC algorithm that provides, at each iteration, a draw of all the trees making up the sum that apporoximates a 
	function f(x). We can call such draw f_d(x), w/ indicating 
	the number of the draw (eg draw 142 out of 500 draw that I want for the target distribution ...?). A prediction is
	usually made by taking the average of all the f_d(x), but a draws from the predictive distribution of a Y given x 
	= x_new can be done by simply picking one d: f_d(x_new) is one draw from the predictive distribution of y given the
	observed xs. Now if we consider Y as a variable under imputation, and X_obs as a set of predictor values corresponding 
	to the observed values on Y (Y_obs), we can obtain a fraw of f_d(x) based on Y_obs and X_obs, and then a draw from
	the posterior predictive distribution of Y_mis given Y_obs, X_obs and the X_mis by picking one f_d(x) and plugging
	in the values of X_mis for each Y_mis.

* **Joint Modelling**
	HeBelin2014 have implemented a a joint modelling multiple imputation algorithm for highdimensional datasets w/ both
	continous and dichotomous variables. They use a generalized multivariate probit model to allow for boht types of variables.
	The idea is that using a generalization of the multivariate probit model could allow to include ordinal as well as nominal
	variables. Citations 10 to 13 in HeBelin2014 could be useful for such purpose. Ideally, include an extension of this 
	paper's proposition that allows for multinomial categorical variables.

** DATA GENERATION **
* For normal data: can follow unstructured covariance matrix generation presented by He Belin 2014
