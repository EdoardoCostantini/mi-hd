# Simulation set up ideas

## Overarching research question
To what extent do the MI imputation methods considered grant statistical validity of the secondary data
analysis when the data imputed is high-dimensional.

## Simulation studies
* 1) multivariate normal data with no structure - Baseline: how well do the methods "recover" infomration on the covaraince structure of a simple unrealistic scenario?
	* factors: 1) % of missing infomration; 2) dimensionality (p<<n, p =~ n, p >> n)
* 2) (multivariate normal) data with latent structure - How do the method perform when missingness is on items that measure latent constructs
	* fixed: % of missing infomration; dimensionality (p >> n);
	* factors: 1) latent construct normality level (+/- deviation from normality: i'm thinking opinions/positions about polirising themes in society); 2) dependency between the latent constructs (+/-, e.g.  if all latent constructs are normally distributed then think of it maybe creating three tiers of correlation levels: low, medium high, where the correaltion coefficients maybe are randomly generated for pairs of constructs between 0 and .33, .33 and .66; .66 and 1)
* 3) mesurement levels - How well do the methods behave with different measurement levels (nominal, ordinal, interval)
	* fixed: % of missing infomration;  dimensionality (p >> n); latent structure
	* factors: 1) percentage of p that are mesaured with each level?
* 4) complex relationships - How well do the methods preserve non-linearity and interactions
	* fixed: % of missing infomration;  dimensionality (p >> n); latent structure?
	* factors: 1) percentage of p that are categorical?
* 5) corssectional EVS datasent with clusters?

## Factorial Design
* fixed: MAR + MCAR; sample size; set of model-relevant, auxiliary and junk variables
* factors: 
	* % missing info; 
	* number of p p<<n, p =~ n, p >> n; 
	* wierd variable distributions (binomials, poissons (e.g. number of children), skewed truncated normals (many 1-10 scales), 
		cut-off (mmm something like age where no one younger than 18 is measured)
	* interaction types and quantity (based maybe on that paper that you read about MICE imputation and interaction types)

## Think about
* y model: there should be (at least) one vairable that is generated according to some model involving 'true' predictors
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

