
# Arbuscular mycorrhizal fungi composition in soils of switchgrass harvested for bioenergy under alternate nitrogen management

Detailed below is the reproducible workflow for the data analysis performed in Dirks & Jackson (2020). All the commands should be executed from the project root directory. 

## 1: Import GitHub Project

Open Terminal and change the current working directory to the location where you want this project to be cloned. Clone the repository with `git clone` and enter the new directory.  

```
git clone https://github.com/aldendirks/amf_metabarcoding.git
cd ./amf_metabarcoding
```

## 2: Acquire Raw Sequence Data

From the project root directory, run the `get_sra.sh` script to install `SRA Toolkit` and download the raw sequence data from NCBI SRA BioProject accession PRJNA590305. After downloading SRA Toolkit, you will be prompted to manually designate a directory to store the sequence files. Use `tab` to navigate to `Cache`, press `o` on the keyboard to select a directory, and use the arrow keys and `enter` to navigate to the `data/seqs/bioinfo/03_css/` directory in the project directory. Make this the SRA download directory.

```
bash scripts/get_sra.sh
```






Alternatively, skip the raw sequence data and just download the processed data files, which total x Mb. 

```
bash scripts/x
```

## 3: Bioinformatics

This third step is only relevant if you downloaded the raw sequence data and want to process them from scratch. Otherwise, skip to step 5, "Analyze AMF Communities".

### 3a: Format Files for `DADA2`

PacBio CCS files are returned in BAM format. To process them with `DADA2` in `R` they first need to be converted to FASTQ format with [`BBMap`](https://sourceforge.net/projects/bbmap/) and the quality scores adjusted to the conventional scale. This script must be run on a Linux OS; I used an Ubuntu virtual machine through VirtualBox.

```
bash scripts/reformat_bam.sh
```

### 3b: Process CCS into ASVs with `DADA2`

Run the `dada.R` script from the command line or open it in R and execute the commands sequentially (recommended). 

```
R scripts/dada2.R
```

### 3c: Cluster ASVs

## 4: Clean Up Data

## 5: Make Phylogenetic Tree

Once all of the data for this project is processed and downloaded, build a phylogenetic tree with the ASVs.

```
```

## 6: Analyze AMF Communities