#!/bin/bash

# Purpose:
# Script to process CCS files with DADA2 using the WEI Scarcity computing cluster.

# From project root directory, usage: 
#   bash scripts/scarcity.sh

####################################################################################################

unset R_LIBS
/opt/bifxapps/R-3.5.3/bin/Rscript scripts/dada2.R