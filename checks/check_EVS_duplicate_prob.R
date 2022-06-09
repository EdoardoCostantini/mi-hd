# Project:   imputeHD-comp
# Objective: Compute probability of a dublicated dataset
# Author:    Anonymized for peer review
# Created:   2022-06-03
# Modified:  2022-06-03

# Initialize resampling study

rm(list=ls())
source("./init_general.R")
source("./exp4_init.R")

# Sample size
N <- 8045
n <- 300

N <- 10
n <- 3


# Sample
ids <- sample(1:N, n, replace = TRUE)

for(r in 1:10){
  for (i in 1:500){
    sample(1:N, n*10, replace = TRUE)
    matrix(sample(1:N, n*10, replace = TRUE), ncol = 10)
  }
}

