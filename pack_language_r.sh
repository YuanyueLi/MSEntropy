#!/bin/bash

# Run R commands to build in bash:
#    devtools::document()

cd language_r
Rscript -e "devtools::document()"
R CMD build .
