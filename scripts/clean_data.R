# Effects of N Treatment on AMF Community Composition in WI Switchgrass Cropping Systems
# Clean Data for Analysis

########################################################################################

# CLEAR ENVIRONMENT
rm(list = ls())

# LOAD LIBRARIES
library(Biostrings); packageVersion("Biostrings") # version 2.52.0
library(picante); packageVersion("picante") # version 1.8
library(plyr); packageVersion("plyr") # version 1.8.4
library(tidyverse); packageVersion("tidyverse") # version 1.2.1
library(vegan); packageVersion("vegan") # version 2.5.5

# SET PATHS
path_cn <- file.path("data", "cn") # path to carbon and nitrogen data
path_glbrc <- file.path("data", "glbrc") # path to GLBRC agronomy and soil data
path_lipids <- file.path("data", "lipids") # path to lipids data
path_master <- file.path("data", "master") # path to master datasets
path_seqs_proc <-  file.path("data", "seqs", "bioinfo", "05_processed") # path to sequence data (to be filtered)
path_seqs_filt <-  file.path("data", "seqs", "bioinfo", "06_filtered") # path to filtered sequence data
path_seqs_sum <- file.path("data", "seqs", "bioinfo", "07_summarized") # path to summarized sequence data
path_tree <- file.path("data", "seqs", "phylogenetics", "raxml", "glomeromycotina") # path to phylogenetic tree

# PROCESS SOIL CARBON NITROGEN DATA
# Import CN data, format, and collapse by experimental unit.
df_cn <- read_csv(file.path(path_cn, "cn_raw.csv"))
df_cn <- df_cn[,-c(3, 5)]
df_cn[,1:3] <- lapply(df_cn[,1:3], as.factor)
df_cn <- df_cn %>% 
  group_by(site, block, nitrogen) %>%
  summarise_all(~mean(.))
df_cn$cn_ratio <- df_cn$per_C / df_cn$per_N
write.csv(df_cn, file = file.path(path_cn, "cn.csv"), row.names = FALSE)

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
write.csv(df_anpp, file = file.path(path_glbrc, "anpp.csv"), row.names=FALSE)
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
write.csv(df_height, file = file.path(path_glbrc, "height.csv"), row.names=FALSE)
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
write.csv(df_soil, file = file.path(path_glbrc, "soil.csv"), row.names=FALSE)
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
write.csv(df_yield, file = file.path(path_glbrc, "yield.csv"), row.names=FALSE)

# PROCESS LIPIDS DATA
# Import lipids data, format, and collapse by experimental unit.
df_lipids <- read_csv(file.path(path_lipids, "lipids_guilds_and_indicators.csv"))
df_lipids <- df_lipids[,-c(3, 15:33)]
df_lipids[,1:3] <- lapply(df_lipids[,1:3], as.factor)
df_lipids <- df_lipids %>% 
  group_by(site, block, nitrogen) %>%
  summarise_all(~mean(.))
df_lipids_amf <- read_csv(file.path(path_lipids, "amf_lipid_ratio.csv"))
df_lipids_amf[,1:3] <- lapply(df_lipids_amf[,1:3], as.factor)
df_lipids_amf <- df_lipids_amf %>% 
  group_by(site, block, nitrogen) %>%
  summarise_all(~mean(.))
write.csv(df_lipids, file = file.path(path_lipids, "lipids.csv"), row.names = FALSE)
write.csv(df_lipids_amf, file = file.path(path_lipids, "lipids_amf.csv"), row.names = FALSE)

# PROCESS ASV MASTER DATA
# Import ASV data
df_asv <- t(read.table(file.path(path_seqs_filt, "ASVs_counts.csv"), header = TRUE, row.names = 1, check.names = FALSE, sep = ","))
df_id <- read.table(file.path(path_seqs_sum, "sample_id.csv"), header = TRUE, row.names = 1, check.names = FALSE, sep = ",")
df_asv <- cbind(df_id[,-3], df_asv)
df_asv <- df_asv %>% 
  group_by(site, block, nitrogen) %>%
  summarise_all(~sum(.))
# Import phylogenetic tree
phy_tree <- read.tree(file.path(path_tree, "RAxML_bipartitions.amf_phylogeny_glomeromycotina"))
rooted_phy_tree <- root(phy_tree, outgroup = 182, resolve.root = TRUE)
rooted_phy_tree_pruned <- prune.sample(df_asv[,4:439], rooted_phy_tree)
phy_dist <- cophenetic(rooted_phy_tree_pruned)
# Make master dataset
df_div <- data.frame(
  "site" = df_asv[,1],
  "block" = df_asv[,2],
  "nitrogen" = df_asv[,3],
  "n" = rowSums(df_asv[,4:439]),
  "rich" = specnumber(df_asv[,4:439]),
  "rich_rar" = rarefy(df_asv[,4:439], 3988),
  "div_shan" = diversity(df_asv[,4:439], index = "shannon"),
  "div_simp" = diversity(df_asv[,4:439], index = "simpson"),
  "even_H" = diversity(df_asv[,4:439], index = "shannon")/log(specnumber(df_asv[,4:439])),
  "faith" = unname(pd(df_asv[,4:439], phy_tree, include.root = FALSE)[1]),
  "mpd" = mpd(as.data.frame(df_asv[,4:439]), phy_dist)
)
df_master <- Reduce(function(...) merge(..., all=TRUE), list(df_anpp[,-5], df_yield[,-5], df_height, df_soil, df_cn, df_lipids, df_lipids_amf, df_div, df_asv))
df_info <- df_master[,1:29]
write.csv(df_asv[,-c(1:3)], file = file.path(path_seqs_sum, "ASVs_counts.csv"), row.names = FALSE)
write.csv(df_info, file = file.path(path_seqs_sum, "sample_info.csv"), row.names = FALSE)
write.csv(df_master, file = file.path(path_master, "master_asv.csv"), row.names = FALSE)

# EXPORT CSV WITH ASV SITE OCCURRENCE
site_occurrence <- df_asv[,-c(2,3)] %>%
  group_by(site) %>%
  summarise_all(~sum(.))
site_occurrence[site_occurrence > 0] <- 1
site_occurrence <- t(site_occurrence)
site_occurrence <- site_occurrence[-1,]
write.csv(site_occurrence, file = file.path(path_seqs_sum, "site_occurrence.csv"), row.names = TRUE)

# PROCESS SPECIES MASTER DATA
tax <- as.matrix(read.table(file.path(path_seqs_filt, "ASVs_tax_custom.csv"), header = TRUE, row.names = 1, check.names = FALSE, sep = ",", na.strings = "unknown"))
df_species <- df_asv[,4:439]
colnames(df_species) <- tax[,8]
df_species <- cbind(sapply(unique(colnames(df_species)[duplicated(colnames(df_species))]), function(x) rowSums(df_species[,grepl(paste(x, "$", sep=""), colnames(df_species))])), df_species[,!duplicated(colnames(df_species)) & !duplicated(colnames(df_species), fromLast = TRUE)])
df_species$site <- df_asv$site; df_species$block <- df_asv$block; df_species$nitrogen <- df_asv$nitrogen 
df_species <- df_species[,c(54:56, 1:53)]
# Make master dataset
df_div <- data.frame(
  "site" = df_species[,1],
  "block" = df_species[,2],
  "nitrogen" = df_species[,3],
  "n" = rowSums(df_species[,4:56]),
  "rich" = specnumber(df_species[,4:56]),
  "rich_rar" = rarefy(df_species[,4:56], 3988),
  "div_shan" = diversity(df_species[,4:56], index = "shannon"),
  "div_simp" = diversity(df_species[,4:56], index = "simpson"),
  "even_H" = diversity(df_species[,4:56], index = "shannon")/log(specnumber(df_species[,4:53]))
)
df_master <- Reduce(function(...) merge(..., all=TRUE), list(df_anpp[,-5], df_yield[,-5], df_height, df_soil, df_cn, df_lipids, df_lipids_amf, df_div, df_species))
write.csv(df_master, file = file.path(path_master, "master_species.csv"), row.names = FALSE)

# PROCESS GENUS MASTER DATA
df_genus <- df_asv[,4:439]
colnames(df_genus) <- tax[,6]
df_genus <- cbind(sapply(unique(colnames(df_genus)[duplicated(colnames(df_genus))]), function(x) rowSums(df_genus[,grepl(paste(x, "$", sep=""), colnames(df_genus))])), df_genus[,!duplicated(colnames(df_genus)) & !duplicated(colnames(df_genus), fromLast = TRUE)])
df_genus$site <- df_asv$site; df_genus$block <- df_asv$block; df_genus$nitrogen <- df_asv$nitrogen 
df_genus <- df_genus[,c(26:28, 1:25)]
# Make master dataset
df_div <- data.frame(
  "site" = df_genus[,1],
  "block" = df_genus[,2],
  "nitrogen" = df_genus[,3],
  "n" = rowSums(df_genus[,4:28]),
  "rich" = specnumber(df_genus[,4:28]),
  "rich_rar" = rarefy(df_genus[,4:28], 3988),
  "div_shan" = diversity(df_genus[,4:28], index = "shannon"),
  "div_simp" = diversity(df_genus[,4:28], index = "simpson"),
  "even_H" = diversity(df_genus[,4:28], index = "shannon")/log(specnumber(df_genus[,4:28]))
)
df_master <- Reduce(function(...) merge(..., all=TRUE), list(df_anpp[,-5], df_yield[,-5], df_height, df_soil, df_cn, df_lipids, df_lipids_amf, df_div, df_genus))
write.csv(df_master, file = file.path(path_master, "master_genus.csv"), row.names = FALSE)

# PROCESS FAMILY MASTER DATA
df_family <- df_asv[,4:439]
colnames(df_family) <- tax[,5]
df_family <- cbind(sapply(unique(colnames(df_family)[duplicated(colnames(df_family))]), function(x) rowSums(df_family[,grepl(paste(x, "$", sep=""), colnames(df_family))])), df_family[,!duplicated(colnames(df_family)) & !duplicated(colnames(df_family), fromLast = TRUE)])
df_family$site <- df_asv$site; df_family$block <- df_asv$block; df_family$nitrogen <- df_asv$nitrogen 
df_family <- df_family[,c(10:12, 1:9)]
# Make master dataset
df_div <- data.frame(
  "site" = df_family[,1],
  "block" = df_family[,2],
  "nitrogen" = df_family[,3],
  "n" = rowSums(df_family[,4:12]),
  "rich" = specnumber(df_family[,4:12]),
  "rich_rar" = rarefy(df_family[,4:12], 3988),
  "div_shan" = diversity(df_family[,4:12], index = "shannon"),
  "div_simp" = diversity(df_family[,4:12], index = "simpson"),
  "even_H" = diversity(df_family[,4:12], index = "shannon")/log(specnumber(df_family[,4:12]))
)
df_master <- Reduce(function(...) merge(..., all=TRUE), list(df_anpp[,-5], df_yield[,-5], df_height, df_soil, df_cn, df_lipids, df_lipids_amf, df_div, df_family))
write.csv(df_master, file = file.path(path_master, "master_family.csv"), row.names = FALSE)

# PROCESS ORDER MASTER DATA
df_order <- df_asv[,4:439]
colnames(df_order) <- tax[,4]
df_order <- cbind(sapply(unique(colnames(df_order)[duplicated(colnames(df_order))]), function(x) rowSums(df_order[,grepl(paste(x, "$", sep=""), colnames(df_order))])), df_order[,!duplicated(colnames(df_order)) & !duplicated(colnames(df_order), fromLast = TRUE)])
df_order$site <- df_asv$site; df_order$block <- df_asv$block; df_order$nitrogen <- df_asv$nitrogen 
df_order <- df_order[,c(5:7, 1:4)]
# Make master dataset
df_div <- data.frame(
  "site" = df_order[,1],
  "block" = df_order[,2],
  "nitrogen" = df_order[,3],
  "n" = rowSums(df_order[,4:7]),
  "rich" = specnumber(df_order[,4:7]),
  "rich_rar" = rarefy(df_order[,4:7], 3988),
  "div_shan" = diversity(df_order[,4:7], index = "shannon"),
  "div_simp" = diversity(df_order[,4:7], index = "simpson"),
  "even_H" = diversity(df_order[,4:7], index = "shannon")/log(specnumber(df_order[,4:7]))
)
df_master <- Reduce(function(...) merge(..., all=TRUE), list(df_anpp[,-5], df_yield[,-5], df_height, df_soil, df_cn, df_lipids, df_lipids_amf, df_div, df_order))
write.csv(df_master, file = file.path(path_master, "master_order.csv"), row.names = FALSE)
