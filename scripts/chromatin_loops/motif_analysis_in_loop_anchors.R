library(GenomicRanges)
library(rtracklayer)
library(universalmotif)
library(Biostrings)
library(monaLisa)
library(ggplot2)
library(ggprism)

# Input data
#############################################
score_quantiles <- c(0.00, 0.25, 0.50, 0.75, 0.80, 0.85, 0.90, 0.95, 0.97, 0.98, 0.99, 1.00)

pka_fn <- "./data/chromatin_loops/data/chromatin_loops/Mlei_loops_annotation_PEU.bed"
gen_fn <- "./data/genome/Mlei_gDNA.fasta"
gix_fn <- "./data/genome/Mlei_gDNA.fasta.fai"
mot_fn <- "Homer_motif_list.motif"
#############################################

# Load top motifs from the de novo HOMER library
#############################################
mot <- homerToPFMatrixList(filename = mot_fn, n = 100L)
pfms <- universalmotif::convert_motifs(mot, class = "TFBSTools-PWMatrix")
mot_m_l <- pfms[[1]]
mot_m_l@ID <- "de_novo_Homer"
mot_m_l@name <- "motif1"

# Prepare genomic bins to scan for the motif
##############################################
gen <- Biostrings::readDNAStringSet(gen_fn, format = "fasta")
gix <- read.table(gix_fn, stringsAsFactors = FALSE)
gix <- gix[,1:2]
colnames(gix) <- c("chr", "length")

gix_seqlengths <- gix$length
names(gix_seqlengths) <- gix$chr
sub_bin <- GenomicRanges::tileGenome(gix_seqlengths, tilewidth = 2000)

sub_gen <- Biostrings::getSeq(gen, sub_bin)
sub_gen <- unlist(sub_gen)
names(sub_gen) <- 1:length(sub_gen)
##############################################

# Load loop regions and genomic background regions
##############################################
pka <- read.csv(pka_fn, header = TRUE, sep = "\t")
pka <- makeGRangesFromDataFrame(pka, keep.extra.columns = TRUE)
pka <- resize(pka, 2000, fix = "center")

seqlevels(pka)  <- gix$chr
seqlengths(pka) <- gix$length
loop_gen <- Biostrings::getSeq(gen, pka)
names(loop_gen) <- pka$name

# Make genomic background bins that exclude loop regions
sub_bin_gr <- unlist(sub_bin)
sub_bin_gr$name <- 1:length(seqnames(sub_bin_gr))
sub_bin_gr <- sub_bin_gr[sub_bin_gr %outside% pka]
bg_gen <- sub_gen[names(sub_gen) %in% sub_bin_gr$name]
##############################################

# Scan genome to identify genome-wide distribution of motif scores
##############################################
sub_gen_r <- monaLisa::findMotifHits(query = mot_m_l,
                                     subject = sub_gen,
                                     min.score = 0,
                                     BPPARAM = BiocParallel::MulticoreParam(4))

# Find max motif score per bin
gw_scores_d <- data.frame(from_bin = sub_gen_r@seqnames, score = sub_gen_r$score)
gw_scores_d_max <- aggregate(score ~ from_bin, data = gw_scores_d, FUN = max)
gw_scores_d_max$code <- "Genome-wide background"

# Find quantile distribution motif scores
options(scipen = 999)
gw_bin_scores_q <- quantile(gw_scores_d_max [ , "score"], score_quantiles)
names(gw_bin_scores_q) <- as.character(score_quantiles)

# Find distribution of scores in loop regions
sub_loop_r <- monaLisa::findMotifHits(query = mot_m_l,
                                      subject = loop_gen,
                                      min.score = 0,
                                      BPPARAM = BiocParallel::MulticoreParam(4))

loop_scores_d <- data.frame(from_bin = sub_loop_r@seqnames, score = sub_loop_r$score)
loop_scores_d_max <- aggregate(score ~ from_bin, data = loop_scores_d, FUN = max)
loop_scores_d_max$code <- "In loop anchors"

# Find distribution of motif scores in background windows
sub_bg_r <- monaLisa::findMotifHits(query = mot_m_l,
                                    subject = bg_gen,
                                    min.score = 0,
                                    BPPARAM = BiocParallel::MulticoreParam(4))
bg_scores_d <- data.frame(from_bin = sub_bg_r@seqnames, score = sub_bg_r$score)
bg_scores_d_max <- aggregate(score ~ from_bin, data = bg_scores_d, FUN = max)
bg_scores_d_max$code <- "Genome-wide background"
##############################################

max_score <- rbind(loop_scores_d_max, bg_scores_d_max)

ggplot(max_score, aes(x = score, fill = code)) +
       geom_density(alpha = 0.7) +
       scale_fill_manual(values = c("whitesmoke", "darkgray")) +
       labs(title = "", x = "Motif score", y = "Fraction of motif sites") +
       scale_x_continuous(n.breaks = 8, guide = "prism_offset", 
                          limits = c(0, 14)) +
       scale_y_continuous(n.breaks = 5, guide = "prism_offset",
                          limits = c(0, 0.50)) +
       geom_vline(xintercept = ( gw_bin_scores_q["0.8"]),
                  colour  = "darkgray", linetype = "dashed") +
       annotate(geom="text", x=(gw_bin_scores_q["0.8"] - 1), y=0.4,
                label=paste("q0.80"), color="darkgray") +
       theme_prism(base_line_size = 0.5, base_family = "sans", base_fontface = "plain") +
       theme(legend.position = "top", legend.text = element_text(margin = margin(r = 20)),
             plot.title = element_text(size = 10), panel.background = element_blank(),
             axis.text.x = element_text(size = 12, hjust = 0.5),
             axis.text.y = element_text(angle = 0, size = 12, hjust = 0.5))
ggsave("GW_Motif_score_EDF9f.pdf", width = 4.7, height = 3.5)



# Scan loop regions for the motif with s motif score above a certain threshold
sub_gen_r <- monaLisa::findMotifHits(query = mot_m_l,
                                     subject = sub_gen,
                                     min.score = gw_bin_scores_q[["0.8"]],
                                     BPPARAM = BiocParallel::MulticoreParam(4))

# Assign genomic coordinates to scaned sequences
sub_gen_df <- as.data.frame(sub_gen_r)
names(sub_gen_df)[1] <- "name"
bindf <- as.data.frame(unlist(sub_bin))
bindf$name <- 1:length(bindf$seqnames)

joindf <- merge(sub_gen_df, bindf, by = "name", all.x = TRUE)
joindf$start <- joindf$start.x + joindf$start.y
joindf$end <- joindf$end.x + joindf$start.y 

joindf <- joindf[,c(10,15,16,5,6,8,9,1)]
colnames(joindf) <- c("seqnames", "start", "end", "motifStrand", "motifSeq", "motifNames", "motifScore", "feature")

# Percentage of loop sites that overlap motif
#############################################
gw080 <- makeGRangesFromDataFrame(joindf, keep.extra.columns = TRUE)

names(pka) <- pka$name
loop_red_ov <- findOverlaps(gw080, pka)
loop_red <- data.frame(gw080[queryHits(loop_red_ov)], pka[subjectHits(loop_red_ov)])
yy <- data.frame("position" = c("in_loop", "outside","total"),
                 "count" = c(length(unique(loop_red$name)),
                 length(unique(pka$name)) - length(unique(loop_red$name)),
                 length(unique(pka$name))))

ggplot(yy[-3,], aes(x = 2, y = count, fill = position)) +
  geom_col(color = "black") +
  coord_polar("y", start = 0) + xlim(c(0.2, 2 + 0.5)) +
  theme_void() +
  labs(title = "Motif presence in loop anchors") +
  geom_text(aes(label = paste(round(count, 1), sep = " ")),
            position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = c("#f8a500eb", "whitesmoke")) +
  theme(legend.text = element_text(size = 12), legend.title = element_blank(),
        text = element_text(size = 12))
ggsave("Loop_anchors_with_motifs_Fig4c.pdf", width = 3, height = 3)