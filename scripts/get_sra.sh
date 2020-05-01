#!/bin/bash

# Purpose:
# Script to download and install SRA toolkit for Mac
# For more information on SRA Toolkit insetallation, or 
#   instructions for other operating systems, see 
#   https://ncbi.github.io/sra-tools/install_config.html

# From the project root directory, usage:
#   bash get_sra.sh

####################################################################################################

# download, unpack, and configure SRA Toolkit
mkdir bin
cd bin
wget "http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-mac64.tar.gz"
tar -xzf sratoolkit.current-mac64.tar.gz
rm sratoolkit.current-mac64.tar.gz
cd ..
export PATH=${PATH}:$PWD/bin/sratoolkit.2.10.5-mac64/bin
mkdir data/seqs/bioinfo
mkdir data/seqs/bioinfo/03_ccs
#vdb-config --interactive

# download files
cd data/seqs/bioinfo/03_ccs
prefetch --type bam --option-file ../../SraAccList.txt
find . -type d -name 'SRR*' -exec bash -c 'mv {}/* "$(dirname {})"' \;
rmdir SRR*