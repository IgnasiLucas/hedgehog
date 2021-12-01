#!/bin/bash

if [ ! -e README.html ]; then
   R --quiet --no-save -e "rmarkdown::render('README.Rmd', output_file = 'README.html')" &> run.log
fi
