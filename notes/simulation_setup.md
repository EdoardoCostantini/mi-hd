# Simulation set up ideas

## Overarching research question
To what extent do the MI imputation methods considered grant statistical validity of the secondary data
analysis when the data imputed is high-dimensional.

## Factorial Design
* fixed: MAR + MCAR; sample size; set of auxiliary and junk variables
* factors: % missing info; number of p p<<n, p \simeq n, p >> n; variable types?

## Think about
* Interactions in the true model: should the types of interactions be factors? (e.g. Doove et al 2014)
	I do not think so as our research question is not related specifically
	to how the different methods handle the presence of interactions. We are
	more focused on which method achieves greater statistical validity of the
	secondary analysis. It could be interesting to simply insert a variety of interaction
	effects based on the results of Doove et al 2014: higher interaction effect size -> higher bias;
	higher correlation between interacting variables -> better statistical performance of estiamtes;
	double ordinal interactions are most well preserved, perfect double disordinal are the hardests.

## Evaluation criteria
* estimates bias
* CI coverage
* FMI (interpreted as the relative loss of efficacy)
* computational speed; complexity of specification 
* complexity of implementation will be good to include both based on your own experience during 
	the project and in general if we find a good comparison meter

