#!/bin/bash

# to run shell script type "./scriptName.sh texFileName.noExt"

# Terminal Workflow for paper 1 manuscript

# Update pdf

# First round of pdflatex
pdflatex $1.tex

# Find all the newly created auxiliary files:
auxFiles=$(find *.aux)

# Run bibtex on each auxiliary file:
for aux in $auxFiles; do bibtex $aux; done

# Final Round of pdflatex
pdflatex $1.tex
pdflatex $1.tex

# Remove useless stuff
rm *.aux
rm *.bbl
rm *.blg
rm *.log
rm *.out

# Moved pdf file to correct location
mv $1.pdf ../pdf/

## NOTES:
## Basic workflow from: http://linuxcommand.org/lc3_wss0010.phpc.w
