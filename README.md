# imputeHD-comp
Repository hosting project high-dimensional imputation comparison.

# Summary
Including a large number of predictors in the imputation model underlying a 
Multiple Imputation (MI) procedure is one of the most challenging tasks 
imputers face.

A variety of high-dimensional MI techniques (HD-MICE) can facilitate this task,
but there has been limited research on their relative performance.
In this study, we investigate a wide range of extant HD-MICE techniques that 
can handle a large number of predictors in the imputation model and general 
missing data patterns.

We assess the relative performance of seven HD-MICE methods with a Monte Carlo
simulation study and a resampling study based on real survey data.
The performance of the methods is defined by the degree to which they facilitate 
unbiased and confidence-valid estimates of the parameters of complete data 
analysis models.

We find that using regularized regression to select the predictors used in the
MI model, and using principal component analysis to reduce the dimensionality 
of auxiliary data produce the best results.

# Contents
This directory contains the following main subfolders:
- code: the main software to run the study
- output: the folder where the results of scripts located in code are stored
- data: where the EVS data should be store after cleaning
- convergence: contains scripts to perform convergence checks
- crossvalidate: contains scripts to perform cross-validation of the ridge penalty for
  one of the methods used in the study (bridge)
- checks: contains scripts checking expected behaviour of different functions
   and set-ups

# How to replicate results

The content of this directory can be used to replicate the results reported in the manuscript: "SMR-21-0138.R1 - High-dimensional imputation for the social sciences: a comparison of state-of-the-art methods"

## Running the simulations

We used R for these simulations.

### Simulation study
1. Installing Dependencies:
   1. Open the script [init_general.R](./code/init_general.R) and install the
      packages with the traditional `install.packages()` function.
   2. Install the package `PcAux` using `devtools::install_github("PcAux")`
   3. Install the package `blasso` by downloading a compatible version of the package from the [package author's website](https://www.asc.ohio-state.edu/hans.11/software/blasso/). If you are running on windows, you need to install g++ to be able to install this package. You can follow these [instructions](https://www3.cs.stonybrook.edu/~alee/g++/g++.html)
   4. Install IVEware by following this [guide](https://www.src.isr.umich.edu/software/iveware/iveware-documentation/installation-guide/)

2. Running the simulation:
   1. Open the script [exp1_init.R](./code/exp1_init.R) and make sure that the parameters and conditions of the simulation study are set to desired values. In particular, pay attention to:
      - `parms$IVEloc` which needs to be set to the correct path for the operating system you are running (for more info look for `~/srclib` [here](https://www.src.isr.umich.edu/software/iveware/iveware-documentation/installation-guide/))
   2. Open the script [exp1_simulation_script_win.R](./code/exp1_simulation_script_win.R)
   3. Make sure the working directory is set to the location of this script (`./code/`)
   4. Define the number of clusters to be used by specifying the first argument in the
      function `makeCluster()`
   5. Run the entire script

### EVS resampling study
1. Installing Dependencies: same as above
2. Preparing the EVS population data: 
   1. Download the EVS 2017 Third pre-release [https://doi.org/10.4232/1.13511](https://doi.org/10.4232/1.13511)
   2. Store it in the data folder inside this project (at the same level as the code 
      folder)
   3. Run the script [exp4_readEVS.R](./code/exp4_readEVS.R) 
3. Running the simulation:
   1. Open the script [exp4_simulation_script_win.R](./code/exp4_simulation_script_win.R)
   2. Make sure the working directory is set to the location of this script (`./code/`)
   3. Define the number of clusters to be used by specifying the first argument in the function 
      `makeCluster()`
   4. Run the entire script

## Obtaining the plots and tables
The procedure is described for the simulation study. 
By using the scripts for "exp4", the same procedure can be followed for the EVS 
resampling study.
1. Open the script [exp1_results.R](./code/exp1_results.R) and make sure you specify 
   the name of the .rds file obtained from the simulation study run.
   This script will extract the results reported in the study.
2. Open the script [exp1_analysis.R](./code/exp1_analysis.R) and make sure you specify 
   the name of the .rds file obtained from the [exp1_results.R](./code/exp1_results.R)
   run.
   To obtain all the plots, you can play around with the parameters defining what 
   is plotted by the script. For example, by changing `pm_grep <- "0.3"` to `0.1`
   you will be able to produce the plots for the smaller proportion of missing 
   cases.

## Keeping track of the results

Because it happens that after getting a review you need to add conditions, repetitions, or tweak other aspects of simulation studies, you need to be able to re-run only certain aspects of the study.
This requires being able to stitch together parts of the results.
Here, I want to keep track of which filenames are important for the results.

### Simulation Study

1. `exp1_simOut_20201130_1006.rds`
   - 1e3 repetitions
   - all the original methods (pre-SMR submission)
2. `exp1_simOut_20220201_1749.rds`
   - 1e3 repetitions
   - only additional methods MI-qp and MI-am run as a result of the SMR review
3. `exp1_simOut_20220201_1749_res.rds`
   - outcome of the exp1_results.R script combining (1) and (2)
4. `exp1_simOut_20220225_1035.rds`
   - 1e3 repetitions
   - re-run of bridge with correct intercept inclusion
5. `exp1_cv_bridge_20220224_1042.rds` 
   - Output for cross-validation of bridge with the correct use of intercept
6. `exp1_simOut_20220225_1035_res.rds`
   - Output for `exp1_results.R` script combining (1), (2), and (4)
7. `exp1_cv_IVEware_20230324_1326.rds`
   - Output for cross-validation of IVEware `minR2` using 70 iterations
8. `exp1_conv_IVEware_20230327_1143.rds`
   - Output for convergence checks for IVEware (above 5 iterations everything seems fine)
   - `exp1_cv_IVEware_20230331_1121.rds` is a version with 70 iterations and 100 multiple imputed datasets
9. `exp1_simOut_20230403_1631.rds`
   - Output for IVEware method
10. `exp1_simOut_20230403_1631_res.rds`
   - Output for `exp1_results.R` script combining (1), (2), (4), and (9)

### Extra Simulation Study on Collinearity

1. `exp1_2_convergence_all_meth_20230403_1027.rds`
- Output for convergence checks for all R native methods.
2. `exp1_2_cv_IVEware_20230405_1715.rds`
- Output for convergence checks for IVEware data.
3. `exp1_2_cv_bridge_20230405_1449.rds`
- Output for cross-validation of `ridge` parameter for bridge
4. `exp1_2_cv_IVEware_20230406_1053.rds`
- Output for cross-validation of `minR2` parameter for IVEware
5. `exp1_2_simOut_20230408_1748.rds`
- 30 repetitions for all methods (contains MI-QP time estimate!)
6. `exp1_2_simOut_20230419_1403.rds`
- 500 for all R-based methods
7. `exp1_2_simOut_20230421_1151.rds`
- 500 repetitions for IVEware method (stepFor / MI-SF)

### Resampling Study
1. `exp4_simOut_20201204_2121.rds`
   - first 500 repetitions
2. `exp4_simOut_20201207_1134.rds`
   - next 500 repetitions
3. `exp4_simOut_20220131_1603.rds`
   - 1e3 repetitions
   - only additional methods MI-qp and MI-am run as a result of the SMR review
4. `exp4_simOut_20220226_0950.rds`
   - 1e3 repetitions
   - re-run of bridge with correct intercept inclusion
5. `exp4_simOut_20230323_1551.rds`
   - 1e3 repetitions
   - run of IVEware with 70 iterations
6. `exp4_simOut_20220226_0950_res.rds`
   - outcome of the exp4_results.R script combining (1), (2), (3), and (4)
7. `exp4_simOut_20230323_1551_res.rds`
   - outcome of the exp4_results.R script combining (1), (2), (3), (4), and (5)
8. `exp4_cv_bridge_20220223_1646.rds`
   - contains the results for cross-validation of bridge with the correct use of intercept
9. `exp4_cv_IVEware_20230322_1841.rds`
   - contains the results for cross-validation of IVEware `minR2` parameter
10. `exp4_cv_IVEware_20230328_1544.rds`
   - contains convergence checks results for IVEware on EVS data