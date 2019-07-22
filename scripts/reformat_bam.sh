#!/bin/bash

# Purpose:
# Script to reformat demultiplexed CCS BAM files to FASTQ using BBMap version 38.50. 

# Information:
# To process CCS with DADA2 in R, files must first be reformated to FASTQ.
# BBMap scripts need to be downloaded from the source and the directory indicated in $PATH. 
# For more information, see [https://sourceforge.net/projects/bbmap/].
# Recommendation: download to /bin, extract tar.gz file, and add to $PATH from
#  project root directory: `export PATH=$PATH:$PWD/bin/bbmap`

# From the project root directory, usage:
#  bash reformat_bam.sh

####################################################################################################

for i in {1..48}
do
   reformat.sh in=data/seqs/bioinfo/03_ccs/ccs_lbc$i.bam out=data/seqs/bioinfo/04_fastq/ccs_lbc$i.fastq
done

for i in {1..9}
do
	mv data/seqs/bioinfo/04_fastq/ccs_lbc$i.fastq data/seqs/bioinfo/04_fastq/ccs_lbc0$i.fastq
done