# Effects of N Treatment on AMF Community Composition in WI Switchgrass Cropping Systems
# Clean Data for Analysis

# Notes:
# Remove the hashtags from the beginning of `write.csv` lines to output cleaned up data

########################################################################################

# CLEAR ENVIRONMENT
rm(list = ls())

# LOAD LIBRARIES
library(Biostrings); packageVersion("Biostrings") # version 2.52.0
library(plyr); packageVersion("plyr") # version 1.8.4
library(tidyverse); packageVersion("tidyverse") # version 1.2.1
library(vegan); packageVersion("vegan") # version 2.5.5

# SET PATHS
path_cn <- file.path("data", "cn") # path to carbon and nitrogen data
path_glbrc <- file.path("data", "glbrc") # path to GLBRC agronomy and soil data
path_lipids <- file.path("data", "lipids") # path to lipids data
path_seqs_proc <-  file.path("data", "seqs", "bioinfo", "05_processed") # path to sequence data (to be filtered)
path_seqs_filt <-  file.path("data", "seqs", "bioinfo", "06_filtered") # path to filtered sequence data
path_seqs_sum <- file.path("data", "seqs", "bioinfo", "07_summarized") # path to summarized sequence data
path_seqs_clust <-  file.path("data", "seqs", "decifr", "uclust", "06_filtered") # path to clustered sequence data
path_master <- file.path("data", "master") # path to master datasets

# PROCESS SOIL CARBON NITROGEN DATA
# Import CN data, format, and collapse by experimental unit.
df_cn <- read_csv(file.path(path_cn, "cn_raw.csv"))
df_cn <- df_cn[,-c(3, 5)]
df_cn[,1:3] <- lapply(df_cn[,1:3], as.factor)
df_cn <- df_cn %>% 
  group_by(site, block, nitrogen) %>%
  summarise_all(~mean(.))
df_cn$cn_ratio <- df_cn$per_C / df_cn$per_N
#write.csv(df_cn, file = file.path(path_cn, "cn.csv"), row.names = FALSE)

# PROCESS GLBRC DATA
# Import aboveground net primary productivity data, format, and collapse by experimental unit. 
df_anpp <- read_csv(file.path(path_glbrc, "588-annual+net+primary+productivity+on+the+marginal+land+experiment.csv"), skip = 31)
df_anpp <- df_anpp[-1,-c(1,3:4,7,8,10:12,15)] 
df_anpp[,1:4] <- lapply(df_anpp[,1:4], as.factor)
df_anpp[,5:6] <- lapply(df_anpp[,5:6], as.numeric)
df_anpp <- df_anpp[grep("Wisconsin", df_anpp$site),] 
df_anpp <- subset(df_anpp, treatment == "G5"); df_anpp <- df_anpp[,-2]
colnames(df_anpp)[2] <- "block"
levels(df_anpp$site) <- c("Lake City", "Escanaba", "Lux Arbor", "Hancock", "Rhinelander", "Oregon")
levels(df_anpp$nitrogen) <- c("control", "amended")
df_anpp$site <- factor(df_anpp$site, levels = c("Hancock", "Oregon", "Rhinelander"))
df_anpp$nitrogen <- factor(df_anpp$nitrogen, levels = c("amended", "control"))
df_anpp <- df_anpp %>%
  group_by(site, block, nitrogen) %>%
  summarise_all(~mean(.))
#write.csv(df_anpp, file = file.path(path_glbrc, "anpp.csv"), row.names=FALSE)
# Import height and lodging data, format, and collapse by experimental unit. 
df_height <- read_csv(file.path(path_glbrc, "592-species+transect+plant+heights.csv"), skip = 25)
df_height <- df_height[-1,-c(5,9)]
df_height[,2:5] <- lapply(df_height[,2:5], as.factor)
df_height$sample_date <- as.Date(df_height$sample_date)
df_height$height <- as.numeric(df_height$height)
df_height <- df_height[grep("G5", df_height$plot),]
df_height <- df_height[df_height$sample_date > "2018-01-01",]; df_height <- df_height[,-1]
df_height <- subset(df_height, site == "Hancock" | site == "Oregon" | site == "Rhinelander")
df_height <- subset(df_height, species == "Panicum virgatum L."); df_height <- df_height[,-4]
colnames(df_height)[2:5] <- c("block", "nitrogen", "height", "lodged")
df_height$block <- gsub("G5-", "", df_height$block)
levels(df_height$nitrogen) <- c("control", "amended")
df_height$site <- factor(df_height$site, levels = c("Hancock", "Oregon", "Rhinelander"))
df_height$nitrogen <- factor(df_height$nitrogen, levels = c("amended", "control"))
df_height <- df_height %>% 
  group_by(site, block, nitrogen) %>%
  summarise_all(~mean(.))
df_height$lodged <- as.numeric(df_height$lodged)
#write.csv(df_height, file = file.path(path_glbrc, "height.csv"), row.names=FALSE)
# Import soil data, format, and collapse by experimental unit. 
df_soil <- read_csv(file.path(path_glbrc, "438-agronomic+soil+chemistry+from+deep+and+surface+soil+cores+marginal+sites.csv"), skip = 36)
df_soil <- df_soil[-1,-c(1, 6:8, 10, 12:15, 18:20)]
df_soil[,c(2:4, 7)] <- lapply(df_soil[,c(2:4, 7)], as.factor)
df_soil[,c(5:6, 8)] <- lapply(df_soil[,c(5:6, 8)], as.numeric)
df_soil <- df_soil[grep("Wisconsin", df_soil$site),]
df_soil <- subset(df_soil, sample_date > "2018-01-01"); df_soil <- df_soil[,-1]
df_soil <- subset(df_soil, treatment == "G5"); df_soil <- df_soil[,-2]
colnames(df_soil)[2] <- "block"
colnames(df_soil)[5] <- "nitrogen"
levels(df_soil$site) <- c("Lake City", "Escanaba", "Lux Arbor", "Hancock", "Rhinelander", "Oregon")
levels(df_soil$nitrogen) <- c("control", "amended")
df_soil$site <- factor(df_soil$site, levels = c("Hancock", "Oregon", "Rhinelander"))
df_soil$nitrogen <- factor(df_soil$nitrogen, levels = c("amended", "control"))
df_soil <- df_soil %>% 
  group_by(site, block, nitrogen) %>%
  summarise_all(~mean(.))
#write.csv(df_soil, file = file.path(path_glbrc, "soil.csv"), row.names=FALSE)
# Import yield data, format, and collapse by experimental unit. 
df_yield <- read_csv(file.path(path_glbrc, "465-dry+matter+yield+marginal+land+experiments.csv"), skip = 29)
df_yield <- df_yield[-1,-c(1,6:7,9,12)]
df_yield[,2:5] <- lapply(df_yield[,2:5], as.factor)
df_yield[,6:7] <- lapply(df_yield[,6:7], as.numeric)
df_yield <- df_yield[grep("Wisconsin", df_yield$site),]
df_yield <- subset(df_yield, harvest_date > "2018-01-01"); df_yield <- df_yield[,-1]
df_yield <- subset(df_yield, treatment == "G5"); df_yield <- df_yield[,-2]
colnames(df_yield)[2:5] <- c("block", "nitrogen", "yield", "moisture")
levels(df_yield$site) <- c("Lake City", "Escanaba", "Lux Arbor", "Hancock", "Rhinelander", "Oregon")
levels(df_yield$nitrogen) <- c("control", "amended")
df_yield$site <- factor(df_yield$site, levels = c("Hancock", "Oregon", "Rhinelander"))
df_yield$nitrogen <- factor(df_yield$nitrogen, levels = c("amended", "control"))
df_yield <- df_yield %>% 
  group_by(site, block, nitrogen) %>%
  summarise_all(~mean(.))
#write.csv(df_yield, file = file.path(path_glbrc, "yield.csv"), row.names=FALSE)

# PROCESS LIPIDS DATA
# Import lipids data, format, and collapse by experimental unit.
df_lipids <- read_csv(file.path(path_lipids, "lipids_guilds_and_indicators.csv"))
df_lipids <- df_lipids[,-c(3, 15:33)]
df_lipids[,1:3] <- lapply(df_lipids[,1:3], as.factor)
df_lipids <- df_lipids %>% 
  group_by(site, block, nitrogen) %>%
  summarise_all(~mean(.))
#write.csv(df_lipids, file = file.path(path_lipids, "lipids.csv"), row.names = FALSE)

# FILTER AMF COMMUNITY COMPOSITION DATA
# Load data
df_asv <- read.table(file.path(path_seqs_proc, "ASVs_counts.tsv"), header = TRUE, row.names = 1, check.names = FALSE, sep = "\t")
df_tax <- read.table(file.path(path_seqs_proc, "ASVs_taxonomy.tsv"), header = TRUE, row.names = 1, check.names = FALSE, sep = "\t")
fas_asv <- readDNAStringSet(file.path(path_seqs_proc, "ASVs.fasta"))
filter <- read.table(file.path(path_seqs_proc, "filter.csv"), header = TRUE, check.names = FALSE, sep = ",")
# Remove non-Glomeromycotina ASVs and ones with incomplete ITS
df_asv <- df_asv[ !(row.names(df_asv) %in% filter[,1]),]
df_asv <- t(df_asv) # transpose matrix so that sites are rows and ASVS are columns
df_tax <- df_tax[ !(row.names(df_tax) %in% filter[,1]),]
fas_asv <- fas_asv[!(names(fas_asv) %in% filter[,1])]
# Write files
#write.csv(df_asv, file = file.path(path_seqs_filt, "ASVs_counts.csv"), row.names = TRUE)
#write.csv(df_tax, file = file.path(path_seqs_filt, "ASVs_taxonomy.csv"), row.names = TRUE)
#writeXStringSet(fas_asv, file = file.path(path_seqs_filt, "ASVs.fasta"), format = "fasta")

# PROCESS ASV MASTER DATA
df_asv <- read.table(file.path(path_seqs_filt, "ASVs_counts.csv"), header = TRUE, row.names = 1, check.names = FALSE, sep = ",")
df_id <- read.table(file.path(path_seqs_filt, "sample_id.csv"), header = TRUE, row.names = 1, check.names = FALSE, sep = ",")
df_asv <- cbind(df_id[,-3], df_asv)
df_asv <- df_asv %>% 
  group_by(site, block, nitrogen) %>%
  summarise_all(~sum(.))
df_div <- data.frame(
  "site" = df_asv[,1],
  "block" = df_asv[,2],
  "nitrogen" = df_asv[,3],
  "div_rich" = specnumber(df_asv[,4:387]),
  "div_rar" = rarefy(df_asv[,4:387], 1949),
  "div_shan" = diversity(df_asv[,4:387], index = "shannon"),
  "div_simp" = diversity(df_asv[,4:387], index = "simpson"))
df_master <- Reduce(function(...) merge(..., all=TRUE), list(df_anpp[,-5], df_yield[,-5], df_height, df_soil, df_cn, df_lipids, df_div, df_asv))
df_info <- df_master[,1:27]
#write.csv(df_asv[,-c(1:3)], file = file.path(path_seqs_sum, "ASVs_counts.csv"), row.names = FALSE)
#write.csv(df_info, file = file.path(path_seqs_sum, "sample_info.csv"), row.names = FALSE)
#write.csv(df_master, file = file.path(path_master, "master_asv.csv"), row.names = FALSE)

# PROCESS OTU MASTER DATA
df_otu_des <- read.table(file.path(path_seqs_clust, "chr1_otus.txt"), header = FALSE, row.names = 1, check.names = FALSE, sep ="", col.names = paste0("V",seq_len(62)), fill = TRUE)
df_otu <- df_asv[,4:387]
for (asv_name in colnames(df_asv[,4:387])){
  colnames(df_otu)[colnames(df_otu) == asv_name] <- names(which(apply(df_otu_des, 1, function(r) any(r %in% asv_name))))
}
df_otu <- cbind(sapply(unique(colnames(df_otu)[duplicated(colnames(df_otu))]), function(x) rowSums(df_otu[,grepl(paste(x, "$", sep=""), colnames(df_otu))])), df_otu[,!duplicated(colnames(df_otu)) & !duplicated(colnames(df_otu), fromLast = TRUE)])
df_div_otu <- data.frame(
  "site" = df_asv[,1],
  "block" = df_asv[,2],
  "nitrogen" = df_asv[,3],
  "div_rich" = specnumber(df_otu),
  "div_rar" = rarefy(df_otu, 1949),
  "div_shan" = diversity(df_otu, index = "shannon"),
  "div_simp" = diversity(df_otu, index = "simpson"))
df_otu <- as.data.frame(cbind(as.matrix(df_asv[,1:3]), as.matrix(df_otu)))
df_master_otu <- Reduce(function(...) merge(..., all=TRUE), list(df_anpp[,-5], df_yield[,-5], df_height, df_soil, df_cn, df_lipids, df_div_otu, df_otu))
colnames(df_master_otu) <- gsub("denovo*", "OTU_", colnames(df_master_otu))
#write.csv(df_master_otu, file = file.path(path_master, "master_otu.csv"), row.names = FALSE)

# PROCESS GENUS-LEVEL MASTER DATA
df_gen <- df_asv[,4:387]
colnames(df_gen) <- df_tax[,6]
df_gen <- df_gen[,!is.na(colnames(df_gen))]
df_gen <- cbind(sapply(unique(colnames(df_gen)[duplicated(colnames(df_gen))]), function(x) rowSums(df_gen[,grepl(paste(x, "$", sep=""), colnames(df_gen))])), df_gen[,!duplicated(colnames(df_gen)) & !duplicated(colnames(df_gen), fromLast = TRUE)])
df_div_gen <- data.frame(
  "site" = df_asv[,1],
  "block" = df_asv[,2],
  "nitrogen" = df_asv[,3],
  "div_rich" = specnumber(df_gen),
  "div_rar" = rarefy(df_gen, 372),
  "div_shan" = diversity(df_gen, index = "shannon"),
  "div_simp" = diversity(df_gen, index = "simpson"))
df_gen <- as.data.frame(cbind(as.matrix(df_asv[,1:3]), as.matrix(df_gen)))
df_master_gen <- Reduce(function(...) merge(..., all=TRUE), list(df_anpp[,-5], df_yield[,-5], df_height, df_soil, df_cn, df_lipids, df_div_gen, df_gen))
#write.csv(df_master_gen, file = file.path(path_master, "master_genus.csv"), row.names = FALSE)

# PROCESS FAMILY-LEVEL MASTER DATA
df_fam <- df_asv[,4:387]
colnames(df_fam) <- df_tax[,5]
df_fam <- df_fam[,!is.na(colnames(df_fam))]
df_fam <- cbind(sapply(unique(colnames(df_fam)[duplicated(colnames(df_fam))]), function(x) rowSums(df_fam[,grepl(paste(x, "$", sep=""), colnames(df_fam))])), df_fam[,!duplicated(colnames(df_fam)) & !duplicated(colnames(df_fam), fromLast = TRUE)])
df_div_fam <- data.frame(
  "site" = df_asv[,1],
  "block" = df_asv[,2],
  "nitrogen" = df_asv[,3],
  "div_rich" = specnumber(df_fam),
  "div_rar" = rarefy(df_fam, 1949),
  "div_shan" = diversity(df_fam, index = "shannon"),
  "div_simp" = diversity(df_fam, index = "simpson"))
df_fam <- as.data.frame(cbind(as.matrix(df_asv[,1:3]), as.matrix(df_fam)))
df_master_fam <- Reduce(function(...) merge(..., all=TRUE), list(df_anpp[,-5], df_yield[,-5], df_height, df_soil, df_cn, df_lipids, df_div_fam, df_fam))
#write.csv(df_master_fam, file = file.path(path_master, "master_family.csv"), row.names = FALSE)

# PROCESS ORDER-LEVEL MASTER DATA
df_ord <- df_asv[,4:387]
colnames(df_ord) <- df_tax[,4]
df_ord <- df_ord[,!is.na(colnames(df_ord))]
df_ord <- cbind(sapply(unique(colnames(df_ord)[duplicated(colnames(df_ord))]), function(x) rowSums(df_ord[,grepl(paste(x, "$", sep=""), colnames(df_ord))])), df_ord[,!duplicated(colnames(df_ord)) & !duplicated(colnames(df_ord), fromLast = TRUE)])
df_div_ord <- data.frame(
  "site" = df_asv[,1],
  "block" = df_asv[,2],
  "nitrogen" = df_asv[,3],
  "div_rich" = specnumber(df_ord),
  "div_rar" = rarefy(df_ord, 1949),
  "div_shan" = diversity(df_ord, index = "shannon"),
  "div_simp" = diversity(df_ord, index = "simpson"))
df_ord <- as.data.frame(cbind(as.matrix(df_asv[,1:3]), as.matrix(df_ord)))
df_master_ord <- Reduce(function(...) merge(..., all=TRUE), list(df_anpp[,-5], df_yield[,-5], df_height, df_soil, df_cn, df_lipids, df_div_ord, df_ord))
#write.csv(df_master_ord, file = file.path(path_master, "master_order.csv"), row.names = FALSE)
