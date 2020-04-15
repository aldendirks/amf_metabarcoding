# Effects of N Treatment on AMF Community Composition in WI Switchgrass Cropping Systems
# Process CCS Files with DADA2

# Notes:
# Before processing with DADA2, PacBio raw sequences must be demultiplexed, assembled into 
#   circular consensus sequences, and converted to FASTQ format
# The following code was adopted from various online tutorials:
#   https://benjjneb.github.io/dada2/tutorial.html
#   https://benjjneb.github.io/LRASManuscript/LRASms_fecal.html
# A discussionon between Alden Dirks and Benjamin Callahan is informative
#   https://github.com/benjjneb/dada2/issues/788

################################################################################################################

# CLEAR ENVIRONMENT
rm(list = ls())

# LOAD LIBRARIES
library(Biostrings); packageVersion("Biostrings") # Version 2.52.0
library(dada2); packageVersion("dada2") # Version 1.14
library(rITSx); packageVersion("rITSx") # Version 0.0.3
library(tidyverse); packageVersion("tidyverse") # Version 1.2.1

# SET PATHS
path_fastq <- file.path("data", "seqs", "bioinfo", "04_fastq")
path_processed <- file.path("data", "seqs", "bioinfo", "05_processed")
path_filtered <- file.path("data", "seqs", "bioinfo", "06_filtered")
path_summarized <- file.path("data", "seqs", "bioinfo", "07_summarized")
path_rds <- file.path("data", "seqs", "rds")
path_ref <- file.path("data", "seqs", "ref")
fns <- list.files(path_fastq, pattern = "fastq", full.names = TRUE)
wSSUmCf <- "TATYGYTCTTNAACGAGGAATC" # amf forward primer
wLSUmBr <- "AACACTCGCAYAYATGYTAGA" # amf reverse primer
rc <- dada2:::rc
print("Paths set")

# READ RDS FILES
# In case you need to pick up where you left off earlier
prim <- readRDS(file.path(path_rds, "prim.rds"))
track <- readRDS(file.path(path_rds, "track.rds"))
drp <- readRDS(file.path(path_rds, "drp.rds"))
err <- readRDS(file.path(path_rds, "err.rds"))
dd <- readRDS(file.path(path_rds, "dd.rds"))
st <- readRDS(file.path(path_rds, "st.rds"))
st_nochim <- readRDS(file.path(path_rds, "st_nochim.rds"))
reten <- readRDS(file.path(path_rds, "reten.rds"))
st_nochim_filt <- readRDS(file.path(path_rds, "st_nochim_filt.rds"))
tax_unite_its <- readRDS(file.path(path_rds, "tax_unite_its.rds"))
tax_rdp_lsu <- readRDS(file.path(path_rds, "tax_rdp_lsu.rds"))
tax_silva_lsu <- readRDS(file.path(path_rds, "tax_silva_lsu.rds"))

# VISUALIZE SEQUENCE QUALITY
plotQualityProfile(fns[1:6])
plotQualityProfile(fns[7:12])
plotQualityProfile(fns[13:18])
plotQualityProfile(fns[19:24])
plotQualityProfile(fns[25:30])
plotQualityProfile(fns[31:36])
plotQualityProfile(fns[37:42])
plotQualityProfile(fns[43:48])

# REMOVE PRIMERS
# Remove primer sequences and orient reverse sequences in the forward direction
# Consider relaxing the number of mismatches allowed in removing primers to max.mismatch = 4 
nops <- file.path(path_fastq, "noprimers", basename(fns))
prim <- removePrimers(fns, nops, primer.fwd = wSSUmCf, primer.rev = dada2:::rc(wLSUmBr), 
                      orient = TRUE, verbose = TRUE, max.mismatch = 2, allow.indels = TRUE)
saveRDS(prim, file.path(path_rds, "prim.rds"))
print("Primers removed")

# LENGTH DISTRIBUTION
# Inspect length distribution of sequences with primers removed to determine how to filter the sequences
lens_fn <- lapply(nops, function(fn) nchar(getSequences(fn)))
lens <- do.call(c, lens_fn)
hist(lens, 100)
mean(lens); median(lens)

# FILTER
# Consider relaxing the minimum quality and maximum expected errors filtering parameters to 
#   minQ = 2 and maxEE = 5
filts <- file.path(path_fastq, "noprimers", "filtered", basename(fns))
track <- filterAndTrim(nops, filts, minQ = 2, minLen = 1300, maxLen = 1800, maxN = 0,
                       rm.phix = FALSE, maxEE = 2, verbose = TRUE)
saveRDS(track, file.path(path_rds, "track.rds"))
print("Sequences filtered")

# LENGTH DISTRIBUTION
# Inspect length distribution of filtered sequences
lens_fn_filts <- lapply(filts, function(fn) nchar(getSequences(fn)))
lens_filts <- do.call(c, lens_fn_filts)
hist(lens_filts, 100)
mean(lens_filts); median(lens_filts)

# DEREPLICATE
# Collapse identical sequences
drp <- derepFastq(filts, verbose = TRUE)
saveRDS(drp, file.path(path_rds, "drp.rds"))
print("Sequences dereplicated")

# LEARN ERRORS
# Errors are unique to each sampling run and must be learned for a sample to properly denoise
err <- learnErrors(drp, errorEstimationFunction = PacBioErrfun, BAND_SIZE = 32, 
                   multithread = TRUE, verbose = TRUE)
plotErrors(err)
saveRDS(err, file.path(path_rds, "err.rds"))
print("Errors learned and relearned")

# DENOISE
# Infer ASV (separate true biological sequences from errors)
dd <- dada(drp, err = err, BAND_SIZE = 32, multithread = TRUE, verbose = TRUE, pool = "pseudo")
saveRDS(dd, file.path(path_rds, "dd.rds"))
print("Sequences denoised")

# SEQUENCE TABLE 
# Generate count table, also called a biome table or ASV matrix
st <- makeSequenceTable(dd); dim(st)
saveRDS(st, file.path(path_rds, "st.rds"))
print("Sequence table created")

# REMOVE CHIMERAS 
# Align sequences with those recovered in greater abundance and determine 
#   which sequences are combinations of others
st_nochim <- removeBimeraDenovo(st, multithread = TRUE, verbose = TRUE)
saveRDS(st_nochim, file.path(path_rds, "st_nochim.rds"))
print("Chimeras removed")

# TRACK RETAINED SEQUENCES
reten <- cbind(ccs = prim[,1], primers = prim[,2], filtered = track[,2], 
               denoised = sapply(dd, function(x) sum(x$denoised)), nonchimeric = apply(st_nochim, 1, sum), 
               retained = apply(st_nochim, 1, sum)/prim[,1])
saveRDS(reten, file.path(path_rds, "reten.rds"))
print("Retained sequences evaluated")

# EXPORT FASTA FILE
# Give our sequence headers more manageable names (ASV_1, ASV_2...) and write out a FASTA of our ASV seqs
asv_seqs <- colnames(st_nochim)
asv_headers <- vector(dim(st_nochim)[2], mode = "character")
for (i in 1:dim(st_nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep = "_")
}
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, file.path(path_processed, "ASVs_processed.fasta"))

# RUN ITSx ON PROCESSED SEQUENCES
# Separate full amplicons into rDNA subunits based on the ITSx algorithm
# ITSx must be installed locally and pointed to in $PATH; see https://microbiology.se/software/itsx/
itsx(in_file = file.path(path_processed, "ASVs_processed.fasta"), out_root = file.path(path_processed, "itsx/itsx"), 
     save_regions = c("all"))

# FILTER PROCCESSED SEQUENCES
# Create multi-locus phylogenetic tree with TBAS in DECIFR
# Create "filter.csv" file to exclude problematic (incomplete rDNA subunits) and non-Glomeromycotina amplicons 
filter <- read.table(file.path(path_filtered, "filter.csv"), header = TRUE, check.names = FALSE, sep = ",")
st_nochim_filt <- st_nochim[,c(!(asv_headers %in% as.character(filter[,1])))]
saveRDS(st_nochim_filt, file.path(path_rds, "st_nochim_filt.rds"))
print("Sequence table (no chimeras) filtered")

# EXPORT FASTA FILE
# Give our sequence headers more manageable names (ASV_1, ASV_2...) and write out a FASTA of our ASV seqs
asv_seqs <- colnames(st_nochim_filt)
asv_headers <- vector(dim(st_nochim_filt)[2], mode = "character")
for (i in 1:dim(st_nochim_filt)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep = "_")
}
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, file.path(path_filtered, "ASVs_filtered.fasta"))

# RUN ITSx ON FILTERED SEQUENCES
# Separate full amplicons into rDNA subunits based on the ITSx algorithm
# ITSx must be installed locally and pointed to in $PATH; see https://microbiology.se/software/itsx/
itsx(in_file = file.path(path_filtered, "ASVs_filtered.fasta"), out_root = file.path(path_filtered, "itsx/itsx"), 
     save_regions = c("all"))

# EXPORT COUNT TABLE
asv_tab <- t(st_nochim_filt)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, file.path(path_filtered, "ASVs_counts.csv"), sep = ",", quote = FALSE, col.names = NA)

# ASSIGN TAXONOMY
# RDP - LSU
st_nochim_lsu <- st_nochim_filt
fas_lsu <- readDNAStringSet(file.path(path_filtered, "itsx/itsx.LSU.fasta"))
colnames(st_nochim_lsu) <- as.character(unname(fas_lsu))
tax_rdp_lsu <- assignTaxonomy(st_nochim_lsu, file.path(path_ref, "RDP_LSU_fixed_train_set_v2.fasta"), 
                              multithread = TRUE, verbose = TRUE, tryRC = TRUE)
saveRDS(tax_rdp_lsu, file.path(path_rds, "tax_rdp_lsu.rds"))
row.names(tax_rdp_lsu) <- sub(">", "", asv_headers)
write.table(tax_rdp_lsu, file.path(path_filtered, "ASVs_tax_rdp_lsu.csv"), sep = ",", 
            quote = FALSE, col.names = NA)
# UNITE - ITS
st_nochim_its <- st_nochim_filt
fas_its <- readDNAStringSet(file.path(path_filtered, "itsx/itsx.full.fasta"))
colnames(st_nochim_its) <- as.character(unname(fas_its))
tax_unite_its <- assignTaxonomy(st_nochim_its, file.path(path_ref, "UNITE_general_release_dynamic_02.02.2019.fasta"), 
                                multithread = TRUE, verbose = TRUE, tryRC = TRUE)
saveRDS(tax_unite_its, file.path(path_rds, "tax_unite_its.rds"))
row.names(tax_unite_its) <- sub(">", "", asv_headers)
write.table(tax_unite_its, file.path(path_filtered, "ASVs_tax_unite_its.csv"), sep = ",", 
            quote = FALSE, col.names = NA)
