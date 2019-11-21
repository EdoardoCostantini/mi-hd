# Simulation set up ideas

## Overarching research question
To what extent do the MI imputation methods considered grant statistical validity of the secondary data
analysis when the data imputed is high-dimensional.

## Factorial Design
* fixed: MAR + MCAR; sample size; set of model-relevant, auxiliary and junk variables
* factors: % missing info; number of p p<<n, p =~ n, p >> n;

## Think about
* y model: there should be (at least) one vairable that is generated according to some model involving 'ture' predictors
	we might work with one linear model and maybe one GLM for categorical vairables. But this depends
	also on whether we want to include only continuous variables in the simulation study or other types 
	as well.
* Interactions in the true model: should we perform a couple of simulation studies, one without and one with 
	interactions (polynomials as well) in the true model? If so, should the types of interactions be 
	factors in the simulation study that includes them? (e.g. Doove et al 2014) I do not think so as our research
	question is not related specifically to how the different methods handle the presence of interactions.
	We are more focused on which method achieves greater statistical validity of the secondary analysis. 
	It could be interesting to simply insert a variety of interaction effects based on the results of Doove
	et al 2014: higher interaction effect size -> higher bias; higher correlation between interacting 
	variables -> better statistical performance of estiamtes; double ordinal interactions are most well
	preserved, perfect double disordinal are the hardests.
* Missingness on auxiliary variables: for now it seems reasonable to assume that missingness occurs in 
	the variables that are true predictors of some y, but not in the auxiliary and junk variables

## Evaluation criteria
* estimates bias
* CI coverage
* FMI (interpreted as the relative loss of efficacy)
* computational speed; complexity of specification

