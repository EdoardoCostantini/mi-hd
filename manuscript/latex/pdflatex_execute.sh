#!/bin/bash

# Terminal Workflow for paper 1 manuscript

# Update pdf

pdflatex 00_manuscript.tex
bibtex 00_manuscript
pdflatex 00_manuscript.tex
pdflatex 00_manuscript.tex

# Remove useless stuff
rm *.aux
rm *.bbl
rm *.blg
rm *.log
rm *.out

# Moved pdf file to correct location
mv 00_manuscript.pdf ../pdf/

## NOTES:
## Basic workflow from: http://linuxcommand.org/lc3_wss0010.phpc.w

