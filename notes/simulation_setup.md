# Simulation set up ideas

## Factorial Design
* fixed: MAR + MCAR; sample size; 
* factors: \% missing info, number of p $p<<n$, $p \simeq n$, $p >> n$; variable types?
* uncertain: completely observed variables/number of variables with missings? Proportions of 
	real auxiliary variables and junk. \textbf{UPDATE} Probably good to have a set of variables
	used in the analytical model that will be imputed and a set of junk/auxiliary variables that
	is large and not afflicted by missing values (you could even predict those with simple 
	regression)

## Evaluation criteria
* estimates bias, 
* CI coverage, 
* FMI (interpreted as the relative loss of efficacy);
* computational speed; complexity of specification 
* complexity of implementation will be good to include both based on your own experience during 
	the project and in general if we find a good comparison meter

