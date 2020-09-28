### Title:    Analysis of results from experiment 2: initialize results for tables
### Project:  Imputing High Dimensional Data
### Author:   Edoardo Costantini
### Created:  2020-07-09
### Notes:    wd should be the code folder

# Session
  rm(list = ls())
  source("./init_general.R")
  library(xtable)

# Read results
  filename <- "exp2_simOut_20200819_1743"
  exp2_res <- readRDS(paste0("../output/", filename, "_res.rds"))
