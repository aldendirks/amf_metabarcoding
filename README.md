
# Arbuscular mycorrhizal fungi composition in soils of switchgrass harvested for bioenergy under alternate nitrogen management

Detailed below is the reproducible workflow for the data analysis performed in Dirks & Jackson (2020). All the commands should be executed from the project root directory. 

## 1: Import GitHub Project

Open Terminal and change the current working directory to the location where you want this project to be cloned. Clone the repository with `git clone` and enter the new directory.  

```
git clone https://github.com/aldendirks/amf_metabarcoding.git
cd ./amf_metabarcoding
```

## 2: Acquire Raw Sequence Data

From the project root directory, run the `get_sra.sh` script to install `SRA Toolkit` and download the raw sequence data in fastq format from NCBI SRA BioProject accession PRJNA590305. You will be prompted to specify a download/cache area for downloaded files. Set the workspace directory to `./data/seqs/bioinfo/04_fastq` by navigating to `Cache`, pressing `o` on the keyboard to select a directory, and using the arrow keys and enter to navigate to the correct location. 

```
bash scripts/get_sra.sh
```

## 3: Bioinformatics

### 3a: Process CCS into ASVs with `DADA2`

Open the `scripts/dada.R` script in `R` and execute the commands sequentially to process the fastq circular consensus sequences into ASVs with `DADA2`. 