library(rtracklayer)
library(GenomicRanges)
library(stringr)
library(ggplot2)
library(ggprism)
library(plyr)
library(reshape2)

####### Functions
#################
tss_tes_annot <- function(gtf_fn, species) {
  gtf <- import(gtf_fn)
  gtf <- gtf[gtf$type == "transcript"]
  gtf$species <- species
  gtf_tss <- resize(gtf, 1, "start")
  gtf_tes <- resize(gtf, 1, "end")
  gtf_tss_df <- data.frame("seqnames" = seqnames(gtf_tss),
                           "start" = start(gtf_tss),
                           "end" = end(gtf_tss),
                           "width" = width(gtf_tss),
                           "strand" = strand(gtf_tss),
                           "gene_id" = gtf_tss$gene_id,
                           "species" = gtf_tss$species)
  gtf_tes_df <- data.frame("seqnames" = seqnames(gtf_tes),
                           "start" = start(gtf_tes),
                           "end" = end(gtf_tes),
                           "width" = width(gtf_tes),
                           "strand" = strand(gtf_tes),
                           "gene_id" = gtf_tes$gene_id,
                           "species" = gtf_tes$species)
  write.table(gtf_tss_df, paste(species, "_TSS_annot.bed", sep = ""),
              col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
  write.table(gtf_tes_df, paste(species, "_TES_annot.bed", sep = ""),
              col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
  gtf_tss$type <- "TSS"
  gtf_tes$type <- "TES"
  gtf <- c(gtf_tss, gtf_tes)
  gtf
}
#################

######### Prepare input data
#########

# AB compartment limits
ab_fn <- list.files("./data/compartmentalization/eigenvalues", "*_region.bed$",
                    full.names = TRUE)

ab <- do.call(rbind, lapply(ab_fn, function(x) {
  df <- read.csv(x, h = FALSE, sep = "\t")
  df$species <- gsub(".*/(.*)\\..*", "\\1", x)
  df$species <- gsub("_.*", "", df$species)
  df
}))

# Identify compartment limits
ab1 <- ab[,c(1,2,4,5)]
ab2 <- ab[,c(1,3,4,5)]
colnames(ab2) <- colnames(ab1)
ab_all <- rbind(ab1, ab2)
ab_all_ag <- aggregate(V4 ~ ., ab_all, FUN = paste, collapse = "")
ab <- ab_all_ag[,c(1,2,2,4,3)]
ab$V2.1 <- ab$V2.1 + 1
ab <- subset(ab, ab$V2 != 0)
ab$V4 <- gsub("IA", "AI", ab$V4)
ab$V4 <- gsub("IB", "BI", ab$V4)
ab$V4 <- gsub("BA", "AB", ab$V4)
ab_gr <- makeGRangesListFromDataFrame(ab,
                                      seqnames.field = "V1",
                                      start.field = "V2",
                                      end.field = "V2.1",
                                      keep.extra.columns = TRUE,
                                      split.field = "species")

# Assign bin sizes used to identify compartments in each species
ab_gr$Cowc$resolution <- 1423
ab_gr$Emue$resolution <- 10790
ab_gr$Mlei$resolution <- 10252
ab_gr$Tadh$resolution <- 4509
ab_gr$Dmel$resolution <- 6878
ab_gr$Hsap$resolution <- 154413


# Loop boundaries
loop_fn <- list.files("./data/chromatin_loops/", "*.bed", full.names = TRUE)

loop <- do.call(rbind, lapply(loop_fn, function(x) {
  df <- read.csv(x, h = FALSE, sep = "\t")
  df$species <- gsub(".*/(.*)\\..*", "\\1", x)
  df$species <- gsub("_.*", "", df$species)
  df
}))

loop_gr <- makeGRangesListFromDataFrame(loop,
                                        seqnames.field = "V1",
                                        start.field = "V2",
                                        end.field = "V3",
                                        keep.extra.columns = TRUE,
                                        split.field = "species")


#TSS and TES gene annotation
gtf_cowc <- tss_tes_annot("./data/gene_annotation_gtf/Cowc_gene_annot.gtf", "Cowc")
gtf_emue <- tss_tes_annot("./data/gene_annotation_gtf/Emue_23chromosomes_gene_annot.gtf", "Emue")
gtf_mlei <- tss_tes_annot("./data/gene_annotation_gtf/Mlei_gene_annot.gtf", "Mlei")
gtf_tadh <- tss_tes_annot("./data/gene_annotation_gtf/Tadh_gene_annot.gtf", "Tadh")
gtf_dmel <- tss_tes_annot("./data/gene_annotation_gtf/Dmel_gene_annot.gtf", "Dmel")
gtf_hsap <- tss_tes_annot("./data/gene_annotation_gtf/Hsap_gene_annot.gtf", "Hsap")

# TSS annotation
tss_fn <- list.files("./", "*_TSS_annot.bed", full.names = TRUE)
tss <- do.call(rbind, lapply(tss_fn, function(x) {
  df <- read.csv(x, h = TRUE, sep = "\t")
  df
}))
tss_gr <- makeGRangesListFromDataFrame(tss,
                                       keep.extra.columns = TRUE,
                                       split.field = "species")

# TES annotation
tes_fn <- list.files("./", "*_TES_annot.bed", full.names = TRUE)
tes <- do.call(rbind, lapply(tes_fn, function(x) {
  df <- read.csv(x, h = TRUE, sep = "\t")
  df
}))
tes_gr <- makeGRangesListFromDataFrame(tes,
                                       keep.extra.columns = TRUE,
                                       split.field = "species")



# H3K4me3 narrow peak
k4me3_fn <- list.files("./data/insulation/H3K4me3", "*.narrowPeak", full.names = TRUE)

k4me3 <- do.call(rbind, lapply(k4me3_fn, function(x) {
  df <- read.csv(x, h = FALSE, sep = "\t")[,1:3]
  df$species <- gsub(".*/(.*)\\..*", "\\1", x)
  df$species <- gsub("_.*", "", df$species)
  df
}))

k4me3_gr <- makeGRangesListFromDataFrame(k4me3,
                                         seqnames.field = "V1",
                                         start.field = "V2",
                                         end.field = "V3",
                                         keep.extra.columns = TRUE,
                                         split.field = "species")


# Insulation boundaries
setwd("./data/insulation/boundary_strength/insulation_score")
filepath <- Sys.glob("*_insulation_score.bg")

insul <- do.call(rbind, lapply(filepath, function(x) {
  df <- read.csv(x, h = FALSE, sep = "\t")
  df$name <- gsub(".*/(.*)\\..*", "\\1", x)
  df
}))

insul$species <- gsub("_.*", "", insul$name)
insul$species <- str_to_title(insul$species)
names(insul)[4] <- "insulation_score"
insul$region <- paste(insul$V1, insul$V2, insul$insulation_score, insul$species, sep = ":")
insul_v2 <- insul
insul_gr <- makeGRangesListFromDataFrame(insul[,-5],
                                         keep.extra.columns = TRUE,
                                         seqnames.field = "V1",
                                         start.field = "V2",
                                         end.field = "V3",
                                         split.field = "species")


### Intersect coordinates of insulation boundaries with other features
#######

# with chromatin loops
IS2loop <- do.call(list, lapply(names(insul_gr), function(i) {
  gr1 <- unlist(insul_gr[i])
  gr2 <- unlist(loop_gr[i])
  gr12_ov <- findOverlaps(gr1, gr2, ignore.strand = TRUE,
                          maxgap = max(unique(width(gr2[i]))))
  gr12 <- unique(gr1[queryHits(gr12_ov)])
  gr12$feature <- paste(names(loop_gr[i]), "Loop", sep = "_")
  gr12
}))
IS2loop <- unlist(as(IS2loop, "GRangesList"))


# with AB compartment limits
IS2ab <- do.call(list, lapply(names(insul_gr), function(i) {
  gr1 <- unlist(insul_gr[i])
  gr2 <- unlist(ab_gr[i])
  gr12_ov <- findOverlaps(gr1, gr2, ignore.strand = TRUE,
                          maxgap = unique(gr2$resolution))
  gr12 <- unique(gr1[queryHits(gr12_ov)])
  gr12$feature <- paste(names(ab_gr[i]), "AB", sep = "_")
  gr12
}))
IS2ab <- unlist(as(IS2ab, "GRangesList"))


# with H3K4me3 peaks
IS2k4me3 <- do.call(list, lapply(names(insul_gr), function(i) {
  gr1 <- unlist(insul_gr[i])
  gr2 <- unlist(k4me3_gr[i])
  gr12_ov <- findOverlaps(gr1, gr2, ignore.strand = TRUE,
                          maxgap = max(unique(width(gr2[i]))))
  gr12 <- unique(gr1[queryHits(gr12_ov)])
  gr12$feature <- paste(names(k4me3_gr[i]), "H3K4me3", sep = "_")
  gr12
}))
IS2k4me3 <- unlist(as(IS2k4me3, "GRangesList"))


# with TSS
IS2tss <- do.call(list, lapply(names(insul_gr), function(i) {
  gr1 <- unlist(insul_gr[i])
  gr2 <- unlist(tss_gr[i])
  gr12_ov <- findOverlaps(gr1, gr2, ignore.strand = TRUE,
                          maxgap = min(unique(width(gr1[i]))) * 2)
  gr12 <- gr1[queryHits(gr12_ov)]
  gr12$feature <- paste(names(tss_gr[i]), "TSS", sep = "_")
  gr12$tss <- gr2$gene_id[subjectHits(gr12_ov)]
  unique(gr12)
}))
IS2tss <- unlist(as(IS2tss, "GRangesList"))


# with TES
IS2tes <- do.call(list, lapply(names(insul_gr), function(i) {
  gr1 <- unlist(insul_gr[i])
  gr2 <- unlist(tes_gr[i])
  gr12_ov <- findOverlaps(gr1, gr2, ignore.strand = TRUE,
                          maxgap = unique(width(gr1[i])))
  gr12 <- gr1[queryHits(gr12_ov)]
  gr12$feature <- paste(names(tes_gr[i]), "TES", sep = "_")
  gr12$tes <- gr2$gene_id[subjectHits(gr12_ov)]
  unique(gr12)
}))
IS2tes <- unlist(as(IS2tes, "GRangesList"))


# Assign features to insulation boundaries
insul <- insul[,6:7]
insul$Loop <- ifelse(insul$region %in% IS2loop$region, 1, 0)
insul$AB <- ifelse(insul$region %in% IS2ab$region, 1, 0)
insul$TSS <- ifelse(insul$region %in% IS2tss$region, 1, 0)
insul$TES <- ifelse(insul$region %in% IS2tes$region, 1, 0)
insul$H3K4me3 <- ifelse(insul$region %in% IS2k4me3$region, 1, 0)
insul$Other <- ifelse(rowSums(insul[,-1:-2]) == 0, 1, 0)
insul_count <- aggregate(. ~ species, insul[,-2], sum)
total <- plyr::count(insul$species)
names(total)[1] <- "species"
insul_count <- merge(insul_count, total, by = "species")
insul_count_perc <- insul_count
insul_count_perc[,-1] <- insul_count_perc[,-1] * 100 / insul_count_perc$freq

# Plot assigned classification as point plot
df <- melt(insul_count_perc[,-8], id.vars = "species")
df$variable <- factor(df$variable, levels = c("Loop", "TSS", "TES", "H3K4me3", "AB", "Other"))
df$species <- factor(df$species, levels = c("Cowc", "Mlei", "Emue", "Tadh", "Dmel", "Hsap"))


ggplot(df, aes(y = species, x = variable, size = value, fill = variable)) +
  geom_point(pch = 21) +
  scale_fill_manual(values = c("#00c1afff","#de99e9ff", "#9b48aaff", "#faca83ff", "#29B6F6", "grey")) +
  labs(x = "", y = "") +
  scale_radius(range = c(2, 12)) +
  theme_minimal(base_line_size = 0.5, base_family = "sans") +
  theme(panel.grid.major = element_line(size = 0.1, colour = "grey"),
        text = element_text(size = 12), legend.title = element_blank(),
        legend.text = element_text(margin = margin(r = 10)),
        axis.text.x = element_text(angle = 0),
        axis.ticks.x = element_blank(), axis.line = element_blank(),
        plot.margin = unit(c(0.8, 0.5, 0.5, 0), "cm"))
ggsave("./Insulation_boundaries_classification_point_plot_EDF5d.pdf", width = 6.5, height = 5)


# Assign features to insulation boundaries using hierarchical classification
insul_v2 <- insul_v2[,6:7]
insul_v2$feature <- ifelse(insul_v2$region %in% IS2loop$region, "Loop",
                           ifelse(insul_v2$region %in% IS2tss$region, "TSS",
                                  ifelse(insul_v2$region %in% IS2ab$region, "AB", "Other")))
insul_v2$count <- 1
insul_v2_ag <- aggregate(count ~ ., insul_v2[,-2], sum)
insul_v2_ag_me$feature <- factor(insul_v2_ag_me$feature,
                                 levels = c("Loop", "TSS", "AB", "Other"))
insul_v2_ag_me$species <- factor(insul_v2_ag_me$species,
                                 levels = c("Cowc", "Mlei", "Emue", "Tadh", "Dmel", "Hsap"))

ggplot(insul_v2_ag_me, aes(x = species, y = count, fill = feature)) +
  geom_bar(position = "fill", stat = "identity", size = 0.1) +
  scale_y_continuous(guide = "prism_offset") +
  scale_fill_manual(values = c("#00c1afff", "#de99e9ff", "#29B6F6", "grey")) +
  labs(x = "", y = "Percentage of annotated IS") +
  theme_prism(base_line_size = 0.5, base_family = "sans", base_fontface = "plain") +
  theme(axis.line.y = element_line(), axis.line.x = element_blank(),
        text = element_text(size = 12),
        legend.title = element_blank(), legend.text = element_text(margin = margin(r = 10)),
        axis.ticks.x = element_blank(), axis.text.x = element_text(angle = 0),
        plot.margin = unit(c(0.8, 0.5, 0.5, 0), "cm"))
ggsave("./Insulation_boundaries_hierarchical_classification_Fig2e.pdf", width = 5.3, height = 3.2)