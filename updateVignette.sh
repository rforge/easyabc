#!/bin/sh

# This script build a new version of the vignette

cd vignettes && \
 svn update && \
 R CMD Sweave EasyABC.Rnw && \
 pdflatex EasyABC.tex && \
 pdflatex EasyABC.tex && \
 cp EasyABC.pdf ../pkg/vignettes/
cd -
