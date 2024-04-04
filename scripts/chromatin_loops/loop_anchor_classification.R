library(rtracklayer)
library(GenomicRanges)
library(ggplot2)
library(ggprism)

setwd("./data/chromatin_loops")

#### To calculate normalized coverage of H3K4me3 and H3K4me2 ChIP signals
#### around center of the loop anchor, estimate required window size based
#### on the width of the ChIP-seq peaks. For M. leidyi the window size is 2 Kb.
################
macs <- read.csv("./data/insulation/H3K4me3/Mlei_H3K4me3_macs_narrow_peaks.narrowPeak",
                 h = FALSE, sep = "\t")
macs$size <- macs$V3 - macs$V2
plot(density(macs$size), xlim = c(0, 5000))

macs2 <- read.csv("./data/insulation/H3K4me2/Mlei_H3K4me2_macs_narrow_peaks.narrowPeak",
                  h = FALSE, sep = "\t")
macs2$size <- macs2$V3 - macs2$V2
plot(density(macs2$size), xlim = c(0, 5000))
################

##### Input data
################
loop_fn <- "./data/chromatin_loops/Mlei_filtered_chromatin_loops_4261.txt"
macsMe2_fn <- "./data/insulation/H3K4me2/Mlei_H3K4me2_macs_narrow_summits.bed"
macsMe3_fn <- "./data/insulation/H3K4me3/Mlei_H3K4me3_macs_narrow_summits.bed"
gw_bin_fn <- "./data/genome/Mlei_gDNA.cooler_makebins_2000bp.bed" # The 2000 bp bins are generated using cooler makebins module.
bedMe2_fn <- "./data/chromatin_loops/H3K4me2_coverage_Mlei_gDNA.cooler_makebins_100bp.bed" # ChIP coverage genome-wide in 100 bp bins.
bedMe3_fn <- "./data/chromatin_loops/H3K4me3_coverage_Mlei_gDNA.cooler_makebins_100bp.bed" # The bins are generated using cooler makebins module.
################


### Input data
##############
# Mark left and right loop anchors, resize all loop anchors to 2000 bp.
loop <- read.csv(loop_fn, h = TRUE, sep = "\t")
l_anchor <- loop[,c(1:3)]
l_anchor$name <- paste("L_loop", seq(1, length(loop$chromosome1)), sep = "_")
r_anchor <- loop[,c(4:6)]
r_anchor$name <- paste("R_loop", seq(1, length(loop$chromosome1)), sep = "_")
colnames(l_anchor) <- colnames(r_anchor)
lr <- rbind(l_anchor, r_anchor)
loop <- makeGRangesFromDataFrame(lr,
                                 seqnames.field = "chromosome1",
                                 start.field = "x1",
                                 end.field = "x2",
                                 keep.extra.columns = TRUE)
loop <- resize(loop, width = 2000, fix = "center")

# Load ChIP-seq and bin data
macsMe2 <- import(macsMe2_fn)
macsMe3 <- import(macsMe3_fn)
gw_bin <- import(gw_bin_fn)

# Load ChIP coverage in 100 bp bins
bedMe2 <- read.csv(bedMe2_fn, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
bedMe2 <- aggregate(V7 ~ ., bedMe2[,-4:-6], max)
names(bedMe2)[4] <- "H3K4me2_raw"

bedMe3 <- read.csv(bedMe3_fn, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
bedMe3 <- aggregate(V7 ~ ., bedMe3[,-4:-6], max)
names(bedMe3)[4] <- "H3K4me3_raw"

chip <- merge(bedMe2, bedMe3, by = c("V1", "V2", "V3"), all = TRUE)
chip$bin <- paste("bin", seq(1, length(chip$V1)), sep = "_")
rownames(chip) <- chip$bin
chip_gr <- makeGRangesFromDataFrame(chip,
                                    keep.extra.columns = TRUE,
                                    seqnames.field = "V1",
                                    start.field = "V2",
                                    end.field = "V3")
##############


# Select background ChIP bins that intersect ChIP summits but not loop anchors
gwGR_bg <- chip_gr[chip_gr %over% c(macsMe2, macsMe3)]
gwGR_bg <- gwGR_bg[gwGR_bg %outside% loop]

# Overlap background bins with 2000 bp bin regions and select the maximum ChIP coverage value.
bg_ov <- findOverlaps(gwGR_bg, gw_bin)
bg <- data.frame(gwGR_bg[queryHits(bg_ov)], gw_bin[subjectHits(bg_ov)], row.names = NULL)
bg_bin <- bg$bin
bg <- bg[,c(9:11, 6, 7, 8)]
bgaMe2 <- aggregate(K4me2_raw ~ ., bg[,c(1:3, 4)], max)
bgaMe3 <- aggregate(K4me3_raw ~ ., bg[,c(1:3, 5)], max)
bga <- merge(bgaMe2, bgaMe3, by = c("seqnames.1", "start.1", "end.1"))
bga$feature <- "Other peaks"
bga$peu <- "control"
bga$loop <- paste("bg", seq(1, length(bga$seqnames.1)), sep = "_")

# Select ChIP coverage bins that intersect loop anchors
tar_ov <- findOverlaps(chip_gr, loop)
tar <- data.frame(chip_gr[queryHits(tar_ov)], loop[subjectHits(tar_ov)], row.names = NULL)
tar <- tar[,c(9:11, 6, 7, 8, 14)]
tarMe2 <- aggregate(K4me2_raw ~ ., tar[,c(1:3, 4, 7)], max)
tarMe3 <- aggregate(K4me3_raw ~ ., tar[,c(1:3, 5, 7)], max)
tarloop <- merge(tarMe2, tarMe3, by = c("seqnames.1", "start.1", "end.1", "name"))
tarloop$feature <- "Loop anchors"
names(tarloop)[4] <- "loop"
tarloop$peu <- ifelse(tarloop$K4me2_raw > tarloop$K4me3_raw, "E", "P")

# Define bins with low ChIP coverage to assign "unknown" feature
plot(density(tarloop$K4me2_raw+ tarloop$K4me3_raw), xlim = c(0,50))
abline(v=15)
tarloop$peu <- ifelse((tarloop$K4me2_raw + tarloop$K4me3_raw) < 15, "U", tarloop$peu)
tarloop$peu <- factor(tarloop$peu, levels = c("P", "E", "U"))

ggplot(tarloop, aes(K4me2_raw, K4me3_raw, color = peu)) +
  geom_point(size = 0.05) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  scale_color_manual(values = c("#F4A460", "#009F6B", "grey80")) +
  scale_x_continuous(limits = c(0, 50), expand = c(0, 0),
                     breaks = seq(0, 50, 10), guide = "prism_offset") +
  scale_y_continuous(limits = c(0, 50), expand = c(0, 0),
                     breaks = seq(0, 50, 10), guide = "prism_offset") +
  labs(title = NULL, x = "H3K4me2, CPM", y = "H3K4me3, CPM") +
  theme_prism(base_line_size = 0.5, base_family = "sans", base_fontface = "plain")  +
  theme(legend.title = element_blank(), legend.position = c(0.8,0.2),  legend.text = element_text(size = 12),
        text = element_text(size = 12),
        axis.line = element_line()) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1)))

# Estimate the number of background ChIP-seq peaks that are not in loop anchors
togGW <- rbind(tarloop, bga)

ggplot(togGW, aes(K4me2_raw, K4me3_raw)) +
  geom_point(aes(colour = feature), size = 0.05, alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 0.4) +
  scale_x_continuous(limits = c(0, 50), expand = c(0, 0),
                     breaks = seq(0, 50, 10), guide = "prism_offset") +
  scale_y_continuous(limits = c(0, 50), expand = c(0, 0),
                     breaks = seq(0, 50, 10), guide = "prism_offset") +
  scale_color_manual(values = c("#00c1afff", "grey80")) +
  labs(title = NULL, x = "H3K4me2, CPM", y = "H3K4me3, CPM") +
  theme_prism(base_line_size = 0.5, base_family = "sans", base_fontface = "plain")  +
  theme(legend.title = element_blank(), legend.position = c(0.8,0.1),  legend.text = element_text(size = 12),
        text = element_text(size = 12),
        axis.line = element_line()) +
  guides(color = guide_legend(override.aes = list(size = 3)))


# Assign enhancer and promoter site annotations to loop anchors
peu_code <- subset(tarloop[,c(4,8)])
names(peu_code)[1] <- "name"
loop <- merge(lr, peu_code, by = "name", all = TRUE)
loop$loop_name <- gsub("L_|R_", "", loop$name)
togl <- subset(loop, grepl("L_.*", loop$name))
togr <- subset(loop, grepl("R_.*", loop$name))

bed12 <- merge(togl[,-1], togr[,-1], by = "loop_name", all = TRUE)
bed12$peu <- paste(bed12$peu.x, bed12$peu.y, sep = "-")
bed12$peu <- ifelse(bed12$peu.x == "U" | bed12$peu.y == "U", "U-U", bed12$peu)
write.table(bed12[,c(2:4, 8:10, 14, 1)], "./data/chromatin_loops/Mlei_loops_annotation_PEU.bed",
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

# Barplot of loop anchor annotations
feature <- bed12$peu
feature <- gsub("E-P", "P-E", feature)
barp <- plyr::count(feature)
barp$percent <- round(barp$freq * 100 / sum(barp$freq), 1)
barp$x <- factor(barp$x, levels = c("P-E", "E-E", "P-P", "U-U"))
barp$species <- "Mlei"

ggplot(barp, aes(fill = x, x = species, y = percent, label = freq)) +
  geom_bar(position = "stack", stat = "identity", size = 0.1) +
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  labs(x = "", y = "Percentage of loops") +
  scale_fill_manual(values = c("#00c1afff", "#5a9cfeff","#faca83ff", "#ccccccff")) +
  scale_y_continuous(expand = c(0,0)) +
  theme_prism(base_line_size = 0.5, base_family = "sans", base_fontface = "plain") +
  theme(axis.line.y = element_line(), axis.line.x = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        text = element_text(size = 12), legend.title = element_blank())
ggsave("./data/chromatin_loops/Mlei_loop_anchor_annotation_Fig2g.pdf", width = 2, height = 3.5)