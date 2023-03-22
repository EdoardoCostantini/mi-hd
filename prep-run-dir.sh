#!/bin/bash
# Project:   mi-spcr
# Objective: Create a run directory for Lisa
# Author:    Edoardo Costantini
# Created:   2022-08-04
# Modified:  2022-08-19

# Define a run name

  runName=$1

# Define a location for the lisa directory

  loc=./runs/$(date +'%Y%m%d')-$runName
  mkdir -p $loc

# Create an empty output folder

  mkdir $loc/output/

# Copy current code folder

  cp -a code $loc/

# Copy the convergence folder

  cp -a convergence $loc/

# Copy the cross-validation folder

  cp -a crossvalidate $loc/

# Copy the input folder

  cp -a input $loc/

# Copy the input folder

  cp -a txt $loc/

# Copy the data folder

  cp -a data $loc/

# Copy the readme

  cp README.md $loc/
