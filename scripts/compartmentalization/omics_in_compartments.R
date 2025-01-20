library("rtracklayer")
library("GenomicRanges")
library("scales")
library("ggprism")
library("ggplot2")

### Functions
#############
percentile_rank <- function(x) round(rank(x) / length(x), 5)

import_omics <- function(file_path, feature) {
  omic_gr <- import(file_path)
  omic_gr$name <- as.numeric(omic_gr$name)
  omic_gr$name <- percentile_rank(omic_gr$name)
  omic_gr$name <- -log2(1 - omic_gr$name)
  omic_gr$omic <- feature
  omic_gr
}

import_rnaseq <- function(file_path, feature) {
  omic <- read.csv(file_path, sep = "\t", h = FALSE)
  names(omic)[5] <- "name"
  omic_gr <- makeGRangesFromDataFrame(omic[,c(1:3,5)],
                                      seqnames.field = "V1",
                                      start.field = "V2",
                                      end.field = "V3",
                                      keep.extra.columns = TRUE)
  omic_gr$name <- as.numeric(omic_gr$name)
  omic_gr$name <- percentile_rank(omic_gr$name)
  omic_gr$name <- -log2(1 - omic_gr$name)
  omic_gr$omic <- feature
  omic_gr
}

region2omics <- function(region, omics) {
  gr_ov <- findOverlaps(region, omics)
  gr <- region[queryHits(gr_ov)]
  gr$score <- as.numeric(omics$name[subjectHits(gr_ov)])
  gr$omics <- omics$omic[subjectHits(gr_ov)]
  df <- data.frame(gr)
  df
}
#############


##### Input AB interval regions
#####
setwd("./data/compartmentalization")
sarc <- import("./eigenvalues/Sarc_compartment_regions_Gaus_distr_region.bed")
cowc <- import("./eigenvalues/Cowc_compartment_regions_Gaus_distr_region.bed")
sros <- import("./eigenvalues/Sros_compartment_regions_Gaus_distr_region.bed")
emue <- import("./eigenvalues/Emue_compartment_regions_Gaus_distr_region.bed")
mlei <- import("./eigenvalues/Mlei_compartment_regions_Gaus_distr_region.bed")
tadh <- import("./eigenvalues/Tadh_compartment_regions_Gaus_distr_region.bed")
nvec <- import("./eigenvalues/Nvec_compartment_regions_Gaus_distr_region.bed")
dmel <- import("./eigenvalues/Dmel_compartment_regions_Gaus_distr_region.bed")
hsap <- import("./eigenvalues/Hsap_compartment_regions_Gaus_distr_region.bed")

regions <- c(sarc, cowc, sros, emue, mlei, tadh, nvec, dmel, hsap)
#####

##### Input omics data for Cowc
#####
k4me3 <- import_omics("./omics/H3K4me3/cowc_chip_H3K4me3_20Kbins.bed", "H3K4me3")
rnaseq <- import_rnaseq("./omics/RNAseq/Cowc_rnaseq_counts_annot.bed", "RNAseq")
te <- read.csv("./omics/TE_density/Cowc_te_density_percentile_rank.bed", h = FALSE, sep = "\t")
gene <- read.csv("./omics/Gene_density/Cowc_gene_density_percentile_rank.bed", h = FALSE, sep = "\t")
colnames(te) <- c("name", "score", "omics")
colnames(gene) <- c("name", "score", "omics")
cowc_omics <- c(k4me3, rnaseq)
cowc_df <- region2omics(cowc, cowc_omics)
cowc_df <- rbind(te, gene, cowc_df[,6:8])
cowc_df$species <- "Cowc"
#####

##### Sarc omics data
#####
k4me3 <- import_omics("./omics/H3K4me3/sarc_chip_H3K4me3_20Kbins.bed", "H3K4me3")
rnaseq <- import_rnaseq("./omics/RNAseq/Sarc_rnaseq_counts_annot.bed", "RNAseq")
te <- read.csv("./omics/TE_density/Sarc_te_density_percentile_rank.bed", h = F, sep = "\t")
gene <- read.csv("./omics/Gene_density/Sarc_gene_density_percentile_rank.bed", h = F, sep = "\t")
colnames(te) <- c("name", "score", "omics")
colnames(gene) <- c("name", "score", "omics")
sarc_omics <- c(k4me3, rnaseq)
sarc_df <- region2omics(sarc, sarc_omics)
sarc_df <- rbind(te, gene, sarc_df[,6:8])
sarc_df$species <- "Sarc"
#####

##### Sros omics data
#####
k4me3 <- import_omics("./omics/H3K4me3/sros_chip_H3K4me3_20Kbins.bed", "H3K4me3")
rnaseq <- import_rnaseq("./omics/RNAseq/Sros_rnaseq_counts_annot.bed", "RNAseq")
te <- read.csv("./omics/TE_density/Sros_te_density_percentile_rank.bed", h = F, sep = "\t")
gene <- read.csv("./omics/Gene_density/Sros_gene_density_percentile_rank.bed", h = F, sep = "\t")
colnames(te) <- c("name", "score", "omics")
colnames(gene) <- c("name", "score", "omics")
sros_omics <- c(k4me3, rnaseq)
sros_df <- region2omics(sros, sros_omics)
sros_df <- rbind(te, gene, sros_df[,6:8])
sros_df$species <- "Sros"
#####

##### Emue omics data
#####
k4me3 <- import_omics("./omics/H3K4me3/emue_chip_H3K4me3_20Kbins.bed", "H3K4me3")
rnaseq <- import_rnaseq("./omics/RNAseq/Emue_rnaseq_counts_annot.bed", "RNAseq")
te <- read.csv("./omics/TE_density/Emue_te_density_percentile_rank.bed", h = FALSE, sep = "\t")
gene <- read.csv("./omics/Gene_density/Emue_gene_density_percentile_rank.bed", h = FALSE, sep = "\t")
colnames(te) <- c("name", "score", "omics")
colnames(gene) <- c("name", "score", "omics")
emue_omics <- c(k4me3, rnaseq)
emue_df <- region2omics(emue, emue_omics)
emue_df <- rbind(te, gene, emue_df[,6:8])
emue_df$species <- "Emue"
#####

##### Mlei omics data
#####
k4me3 <- import_omics("./omics/H3K4me3/mlei_chip_H3K4me3_20Kbins.bed", "H3K4me3")
rnaseq <- import_rnaseq("./omics/RNAseq/Mlei_rnaseq_counts_annot.bed", "RNAseq")
te <- read.csv("./omics/TE_density/Mlei_te_density_percentile_rank.bed", h = FALSE, sep = "\t")
gene <- read.csv("./omics/Gene_density/Mlei_gene_density_percentile_rank.bed", h = FALSE, sep = "\t")
colnames(te) <- c("name", "score", "omics")
colnames(gene) <- c("name", "score", "omics")
mlei_omics <- c(k4me3, rnaseq)
mlei_df <- region2omics(mlei, mlei_omics)
mlei_df <- rbind(te, gene, mlei_df[,6:8])
mlei_df$species <- "Mlei"
#####

##### Tadh omics data
#####
k4me3 <- import_omics("./omics/H3K4me3/tadh_chip_H3K4me3_20Kbins.bed", "H3K4me3")
rnaseq <- import_rnaseq("./omics/RNAseq/Tadh_rnaseq_counts_annot.bed", "RNAseq")
te <- read.csv("./omics/TE_density/Tadh_te_density_percentile_rank.bed", h = FALSE, sep = "\t")
gene <- read.csv("./omics/Gene_density/Tadh_gene_density_percentile_rank.bed", h = FALSE, sep = "\t")
colnames(te) <- c("name", "score", "omics")
colnames(gene) <- c("name", "score", "omics")
tadh_omics <- c(k4me3, rnaseq)
tadh_df <- region2omics(tadh, tadh_omics)
tadh_df <- rbind(te, gene, tadh_df[,6:8])
tadh_df$species <- "Tadh"
#####

##### Nvec omics data
#####
k4me3 <- import_omics("./omics/H3K4me3/nvecElav_chip_H3K4me3_20Kbins.bed", "H3K4me3")
rnaseq <- import_rnaseq("./omics/RNAseq/Nvec_rnaseq_counts_annot.bed", "RNAseq")
te <- read.csv("./omics/TE_density/Nvec_te_density_percentile_rank.bed", h = F, sep = "\t")
gene <- read.csv("./omics/Gene_density/Nvec_gene_density_percentile_rank.bed", h = F, sep = "\t")
colnames(te) <- c("name", "score", "omics")
colnames(gene) <- c("name", "score", "omics")

nvec_omics <- c(k4me3, rnaseq)
nvec_df <- region2omics(nvec, nvec_omics)
nvec_df <- rbind(te, gene, nvec_df[,6:8])
nvec_df$species <- "Nvec"
#####

##### Dmel omics data
#####
k4me3 <- import_omics("./omics/H3K4me3/dmel_chip_H3K4me3_20Kbins.bed", "H3K4me3")
rnaseq <- import_rnaseq("./omics/RNAseq/Dmel_rnaseq_counts_annot.bed", "RNAseq")
te <- read.csv("./omics/TE_density/Dmel_te_density_percentile_rank.bed", h = FALSE, sep = "\t")
gene <- read.csv("./omics/Gene_density/Dmel_gene_density_percentile_rank.bed", h = FALSE, sep = "\t")
colnames(te) <- c("name", "score", "omics")
colnames(gene) <- c("name", "score", "omics")
dmel_omics <- c(k4me3, rnaseq)
dmel_df <- region2omics(dmel, dmel_omics)
dmel_df <- rbind(te, gene, dmel_df[,6:8])
dmel_df$species <- "Dmel"
#####

##### Hsap omics data
#####
k4me3 <- import_omics("./omics/H3K4me3/hsap_chip_H3K4me3_20Kbins.bed", "H3K4me3")
rnaseq <- import_rnaseq("./omics/RNAseq/Hsap_rnaseq_counts_annot.bed", "RNAseq")
te <- read.csv("./omics/TE_density/Hsap_te_density_UCSC_percentile_rank.bed", h = FALSE, sep = "\t")
gene <- read.csv("./omics/Gene_density/Hsap_gene_density_percentile_rank.bed", h = FALSE, sep = "\t")
colnames(te) <- c("name", "score", "omics")
colnames(gene) <- c("name", "score", "omics")
hsap_omics <- c(k4me3, rnaseq)
hsap_df <- region2omics(hsap, hsap_omics)
hsap_df <- rbind(te, gene, hsap_df[,6:8])
hsap_df$species <- "Hsap"
#####

df <- rbind(sarc_df, cowc_df, sros_df, emue_df, mlei_df, nvec_df, tadh_df, dmel_df, hsap_df)

df$species <- factor(df$species,
                     levels = c("Sarc", "Cowc", "Sros", "Mlei", "Emue", "Tadh", "Nvec", "Dmel", "Hsap"))
df$omics <- factor(df$omics,
                  levels = c("H3K4me3", "RNAseq", "TE_density", "Gene_density"),
                  labels = c("H3K4me3", "RNAseq", "TE density", "Gene density"))
df$name <- factor(df$name,
                  levels = c("A", "I", "B"))

ggplot(df, aes(x = species, y = score, fill = name)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = c("#f56652", "#fed976",  "#0060d2")) +
  labs(x = "", y = "-log2(1 - signal quantile)") +
  scale_y_continuous(limits = c(0, 7), guide = "prism_offset") +
  theme_prism(base_line_size = 0.5, base_family = "sans", base_fontface = "plain")  +
  theme(legend.position = "top",
        text = element_text(size = 12),
        axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
  facet_grid(omics ~ species, scales = "free", space = "free_x", switch = "x") +
  theme(strip.text = element_text(size = 12), strip.background.y = element_rect())
ggsave("Compartments_vs_Omics_EDF3d.pdf", width = 8, height = 4.5)