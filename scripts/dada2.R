# Effects of N Treatment on AMF Community Composition in WI Switchgrass Cropping Systems
# Process CCS Files with DADA2

# Notes:
# Before processing with DADA2, sequences must be demultiplexed and converted to FASTQ format
# The following code was adopted from various online tutorials:
#   https://benjjneb.github.io/dada2/tutorial.html
#   https://benjjneb.github.io/LRASManuscript/LRASms_fecal.html
# A discussionon between Alden Dirks and Benjamin Callahan is informative
#   https://github.com/benjjneb/dada2/issues/788

########################################################################################

# CLEAR ENVIRONMENT
rm(list = ls())

# LOAD LIBRARIES
library(dada2); packageVersion("dada2") # Version 1.12.1

# SET PATHS
path_fastq <- file.path("data", "seqs", "bioinfo", "04_fastq")
path_processed <- file.path("data", "seqs", "bioinfo", "05_processed")
path_rds <- file.path("data", "seqs", "rds")
path_ref <- file.path("data", "seqs", "ref")
fns <- list.files(path_fastq, pattern = "fastq", full.names = TRUE)
wSSUmCf <- "TATYGYTCTTNAACGAGGAATC" # amf forward primer
wLSUmBr <- "AACACTCGCAYAYATGYTAGA" # amf reverse primer
rc <- dada2:::rc
print("Paths set")

# READ RDS FILES
# In case you need to pick up where you left off earlier
#prim <- readRDS(file.path(path_rds_amf, "prim.rds"))
#track <- readRDS(file.path(path_rds_amf, "track.rds"))
#drp <- readRDS(file.path(path_rds_amf, "drp.rds"))
#err <- readRDS(file.path(path_rds_amf, "err.rds"))
#dd <- readRDS(file.path(path_rds_amf, "dd.rds"))
#st <- readRDS(file.path(path_rds_amf, "st.rds"))
#st_nochim <- readRDS(file.path(path_rds_amf, "st_nochim.rds"))
#tax <- readRDS(file.path(path_rds, "tax.rds"))
#reten <- readRDS(file.path(path_rds, "reten.rds"))

# VISUALIZE SEQUENCE QUALITY
# Run if using in R
#plotQualityProfile(fns[1:6])
#plotQualityProfile(fns[7:12])
#plotQualityProfile(fns[13:18])
#plotQualityProfile(fns[19:24])
#plotQualityProfile(fns[25:30])
#plotQualityProfile(fns[31:36])
#plotQualityProfile(fns[37:42])
#plotQualityProfile(fns[43:48])

# REMOVE PRIMERS
# Remove primer sequences and orient reverse sequences in the forward direction
# Consider relaxing the number of mismatches allowed in removing primers: max.mismatch = 4 
nops <- file.path(path_fastq, "noprimers", basename(fns))
prim <- removePrimers(fns, nops, primer.fwd = wSSUmCf, primer.rev = dada2:::rc(wLSUmBr), orient = TRUE, verbose = TRUE, max.mismatch = 2)
saveRDS(prim, file.path(path_rds, "prim.rds"))
print("Primers removed")

# LENGTH DISTRIBUTION
# Inspect length distribution of sequences with primers removed to determine how to filter the sequences
lens_fn <- lapply(nops, function(fn) nchar(getSequences(fn)))
lens <- do.call(c, lens_fn)
hist(lens, 100)
mean(lens); median(lens)

# FILTER
# Consider relaxing the minimum quality and maximum expected errors filtering parameters: minQ = 2, maxEE = 5
filts <- file.path(path_fastq, "noprimers", "filtered", basename(fns))
track <- filterAndTrim(nops, filts, minQ = 3, minLen = 1000, maxLen = 1600, maxN = 0, rm.phix = FALSE, maxEE = 2, verbose = TRUE)
saveRDS(track, file.path(path_rds, "track.rds"))
print("Sequences filtered")

# DEREPLICATE
# Collapse identical sequences
drp <- derepFastq(filts, verbose = TRUE)
saveRDS(drp, file.path(path_rds, "drp.rds"))
print("Sequences dereplicated")

# LEARN ERRORS
# Errors are unique to each sampling run and must be learned for a sample to properly denoise
err <- learnErrors(drp, errorEstimationFunction = PacBioErrfun, BAND_SIZE = 32, multithread = TRUE, verbose = TRUE)
plotErrors(err)
saveRDS(err, file.path(path_rds, "err.rds"))
print("Errors learned and relearned")

# DENOISE
# Infer ASVs (separate true biological sequences from errors)
dd <- dada(drp, err = err, BAND_SIZE = 32, multithread = TRUE, verbose = TRUE, pool = "pseudo")
saveRDS(dd, file.path(path_rds, "dd.rds"))
print("Sequences denoised")

# SEQUENCE TABLE 
# Generate count table, also called a biome table or ASV matrix
st <- makeSequenceTable(dd); dim(st)
saveRDS(st, file.path(path_rds, "st.rds"))
print("Sequence table created")

# REMOVE CHIMERAS 
# Align sequences with those recovered in greater abundance and determine which sequences are combinations of others
st_nochim <- removeBimeraDenovo(st, multithread = TRUE, verbose = TRUE)
saveRDS(st_nochim, file.path(path_rds, "st_nochim.rds"))
print("Chimeras removed")

# TRACK RETAINED SEQUENCES
reten <- cbind(ccs = prim[,1], primers = prim[,2], filtered = track[,2], denoised = sapply(dd, function(x) sum(x$denoised)), nonchimeric = apply(st_nochim, 1, sum), retained = apply(st_nochim, 1, sum)/prim[,1])
saveRDS(reten, file.path(path_rds, "reten.rds"))
print("Retained sequences evaluated")

# ASSIGN TAXONOMY
tax <- assignTaxonomy(st_nochim, file.path(path_ref, "custom_ref_dada2.fasta"), multithread = TRUE, verbose = TRUE, tryRC = TRUE)
saveRDS(tax, file.path(path_rds, "tax.rds"))

# EXPORT FILES
# Give our sequence headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(st_nochim)
asv_headers <- vector(dim(st_nochim)[2], mode = "character")
for (i in 1:dim(st_nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep = "_")
}
# Making and writing out a FASTA of our final ASV seqs
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, file.path(path_processed, "ASVs.fasta"))
# Create count table
asv_tab <- t(st_nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, file.path(path_processed, "ASVs_counts.tsv"), sep = "\t", quote = FALSE, col.names = NA)
# Crate taxonomy table
row.names(tax) <- sub(">", "", asv_headers)
write.table(tax, file.path(path_processed, "ASVs_taxonomy.tsv"), sep = "\t", quote = FALSE, col.names = NA)
