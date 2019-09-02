# Effects of N Treatment on AMF Community Composition in WI Switchgrass Cropping Systems

Detailed below is the reproducible workflow for the data analysis performed in (to be published). Note: all the commands below should be executed from the project root directory. 

## 1: Import GitHub Project

Open Terminal and change the current working directory to the location where you want this project to be cloned. Clone the repository with `git clone`.  

```
git clone https://github.com/aldendirks/amf_metabarcoding.git
```

## 2: Acquire Data

From the project root directory, run the `x` script to download all data (raw and processed) totalling x Gb.

```
bash scripts/x
```

Alternatively, skip the raw sequence data and just download the processed data files, which total x Mb. 

```
bash scripts/x
```

## 3: Bioinformatics

This third step is only relevant if you downloaded the raw sequence data and want to process them from scratch. Otherwise, skip to step 5, "Analyze AMF Communities".

### 3a: Format Files for `DADA2`

PacBio CCS files are returned in BAM format. To process them with `DADA2` in `R` they first need to be converted to FASTQ format with [`BBMap`](https://sourceforge.net/projects/bbmap/) and the quality scores adjusted to the conventional scale. This script must be run on a Linux OS; I used an Ubuntu virutal machine through VirtualBox.

```
bash scripts/reformat_bam.sh
```

### 3b: Process CCS into ASVs with `DADA2`

I processed the CCS files on a 256 GB RAM computing cluster through the Wisconsin Energy Institute (called Scarcity). The script below is just a call to `dada2.R`, unsetting some Scarcity-specific parameters to run properly.

```
bash scripts/scarcity.sh
```

Alternatively, run `dada2.R` directly in the command line or open it in `R` and follow along. 

### 3c: Cluster ASVs

## 4: Clean Up Data

## 5: Analyze AMF Communities