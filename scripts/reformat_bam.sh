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

# download and unpack BBMap
cd bin
wget "https://sourceforge.net/projects/bbmap/files/BBMap_38.82.tar.gz"
tar -xzf BBMap_38.82.tar.gz
rm BBMap_38.82.tar.gz
export PATH=${PATH}:$PWD/bbmap

# download and unpack SAM tools (to process BAM files)
wget "https://sourceforge.net/projects/samtools/files/samtools-1.10.tar.bz2"
tar -xzf samtools-1.10.tar.bz2
rm samtools-1.10.tar.bz2
export PATH=${PATH}:$PWD/samtools-1.10
cd ..

# reformat BAM files 
mkdir data/seqs/bioinfo/fastq
for i in {1..9}
do
   reformat.sh in=data/seqs/bioinfo/bam/ccs_lbc0$i.bam out=data/seqs/bioinfo/fastq/ccs_lbc0$i.fastq
done
for i in {10..48}
do
   reformat.sh in=data/seqs/bioinfo/bam/ccs_lbc$i.bam out=data/seqs/bioinfo/fastq/ccs_lbc$i.fastq
done

reformat.sh in=data/seqs/bioinfo/bam/ccs_lbc01.bam out=data/seqs/bioinfo/fastq/ccs_lbc01.fastq