#!/bin/bash

## Project:	parallel_example
## Title:	Create Stopos Pool
## Author:	Edoardo Costantini
## Created:	2021-02-03
## Notes:	Assumes you are in the code folder of the project

## USAGE on LISA:
##   ./pool_script PATH_TO_POOL_LINES
##
## ARGS:
##   PATH_TO_POOL_LINES = Path to the text file with the number of lines to use for
##			  the stopos pool

# Create Stopos Pool
stopos create -p pool	# to create an empty pool of parameters
stopos -p pool add "$1"	# to put the parameters as lines in the pool
stopos status		# print a description of the resulting pool

