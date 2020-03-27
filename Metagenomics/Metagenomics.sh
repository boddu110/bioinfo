#!/usr/bin/env bash
# renderMerge.sh
R -e "rmarkdown::render('Metagenomics.Rmd', output_format='all')"
