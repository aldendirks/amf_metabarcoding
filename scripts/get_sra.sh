#!/bin/bash

# Purpose:
# Script to download and install SRA toolkit for Mac and download raw PacBio circular consensus 
#   sequences BAM files
# For more information on SRA Toolkit installation, or instructions for other operating systems, see 
#   https://ncbi.github.io/sra-tools/install_config.html

# From the project root directory, usage:
#   bash get_sra.sh

####################################################################################################

# download and unpack SRA Toolkit
cd bin
wget "http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-mac64.tar.gz"
tar -xzf sratoolkit.current-mac64.tar.gz
rm sratoolkit.current-mac64.tar.gz
export PATH=${PATH}:$PWD/sratoolkit.2.10.5-mac64/bin
vdb-config -i
cd ..

# download files
cd data/seqs/bioinfo/04_fastq
prefetch --option-file ../../ncbi/SraAccList.txt
cd sra
fasterq-dump --split-files *
cd ..
mv sra/*.fastq .
rm -r sra

# rename files
n=2
for file in *.sra.fastq
do
  mv $file $(sed -n ${n}p ../../ncbi/SraRunTable.txt | cut -d ',' -f 37).fastq
  n=$((n+1))
done