#!/bin/bash

#SBATCH --job-name=raxml_acd
#SBATCH --mail-user=adirks@umich.edu
#SBATCH --mail-type=BEGIN,END
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu=2000m 
#SBATCH --time=120:00:00
#SBATCH --account=tyjames
#SBATCH --partition=standard
#SBATCH --output=/home/%u/%x-%j.log

# Load modules
module load Bioinformatics
module load raxml

# The application(s) to execute along with its input arguments and options:
#   raxmlHPC: command for sequential and parallel Maximum Likelihood based inference of large phylogenetic trees
#   -s: sequence file name (aligned FASTA or phylip format)
#   -n: output file name (not path, prefix for output files)
#   -m: substitution model
#   -p: parsimony random seed (don't use 12345, as in the examples!)
raxmlHPC-PTHREADS -s ~/projects/amf_metabarcoding/data/seqs/phylogenetics/mafft/amf_phylogeny_mafft.fasta -n amf_phylogeny -m GTRGAMMA -p 19930327 -f a -# 100 -T 8 -x 19930327