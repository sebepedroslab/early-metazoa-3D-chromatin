library("GenomicRanges")
library("ggplot2")
library("ggprism")

### Functions
#############
overlap_IS2window <- function(window_gr, IS_gr) {
  gr_ov <- findOverlaps(window_gr, IS_gr)
  IS <- window_gr[queryHits(gr_ov)]
  IS$score <- IS_gr$score[subjectHits(gr_ov)]
  IS$window <- IS_gr$V5[subjectHits(gr_ov)]
  ISdf <- data.frame(IS)
  ISdf
}

merge_w1_w2 <- function(w1, w2, name) {
  ww <- rbind(w1, w2)
  threshold <- quantile(ww$boundary_strength, probs = 0.25)
  
  ww <- subset(ww, ww$boundary_strength > threshold)
  interval_gr <- reduce(makeGRangesFromDataFrame(ww,
                                          seqnames.field = "chrom",
                                          start.field = "start",
                                          end.field = "end"))
  ww_gr <- makeGRangesFromDataFrame(ww, keep.extra.columns = TRUE,
                                    seqnames.field = "chrom",
                                    start.field = "start",
                                    end.field = "end")
  ovgr <- findOverlaps(interval_gr, ww_gr)
  ovdf <- data.frame(interval_gr[queryHits(ovgr)], ww_gr[subjectHits(ovgr)])
  IS <- aggregate(log2_IS ~ ., ovdf[,c(1:3,11)], min)
  boundaries <- aggregate(boundary_strength ~ ., ovdf[,c(1:3,12)], max)
  write.table(IS, paste("./insulation_score/", name, "_insulation_score.bg", sep = ""),
              col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
  write.table(boundaries, paste("./boundary_strength/", name, "_boundary_strength.bg", sep = ""),
              col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
  df <- cbind(IS, boundaries)
  df
}
#############


##### Input insulation score values (weak and strong boundaries)
#####
setwd("./data/insulation")

ISfn <- list.files("./mlei_bg", "mlei.*", full.names = TRUE)
IS <- do.call(rbind, lapply(ISfn, function(x) {
  df <- read.csv(x, h = FALSE, sep = "\t")
  df$name <- gsub(".*/(.*)\\..*", "\\1", x)
  df$V4 <- as.numeric(gsub("NaN", NA, df$V4))
  df <- na.omit(df)
  names(df)[4] <- "score"
  df
}))

IS_gr <- makeGRangesListFromDataFrame(IS, split.field = "name",
                                      keep.extra.columns = TRUE,
                                      seqnames.field = "V1",
                                      start.field = "V2",
                                      end.field = "V3")
#####

### Input +/- 15 kb window regions around strong insulation boundaries
#####
windowfn <- list.files("./metaplot", "mlei.*", full.names = TRUE)

window <- do.call(rbind, lapply(windowfn, function(x) {
  df <- read.csv(x, h = FALSE, sep = "\t")
  df$res <- gsub(".*/(.*)\\..*", "\\1", x)
  df$res <- gsub("_strong_boundary_region_slop_windows", "", df$res)
  df
}))

window_gr <- makeGRangesListFromDataFrame(window, split.field = "res",
                                          keep.extra.columns = TRUE,
                                          seqnames.field = "V1",
                                          start.field = "V2",
                                          end.field = "V3")
#####

### Calculate IS around boundaries with a window of +/- 15 kb
is1 <- overlap_IS2window(window_gr[[1]], IS_gr[[1]])
is1$name <- names(window_gr[1])

is2 <- overlap_IS2window(window_gr[[2]], IS_gr[[2]])
is2$name <- names(window_gr[2])

is3 <- overlap_IS2window(window_gr[[3]], IS_gr[[3]])
is3$name <- names(window_gr[3])

is4 <- overlap_IS2window(window_gr[[4]], IS_gr[[4]])
is4$name <- names(window_gr[4])

is1234 <- rbind(is1, is2, is3, is4)

### Prepare for plotting
ISdf_score <- aggregate(score ~ ., is1234[,6:9], mean)
ISdf_score$V4 <- as.numeric(ISdf_score$V4)

ISdf_score$name <- factor(ISdf_score$name,
                          levels = c("mlei_400bpRes", "mlei_1000bpRes", "mlei_2000bpRes", "mlei_4000bpRes"),
                          labels = c("bin 400bp", "bin 1000bp", "bin 2000bp", "bin 4000bp"))

ISdf_score$window <- factor(ISdf_score$window, levels = c("wx5", "wx10", "wx25"),
                            labels = c("window of 5*bin", "window of 10*bin", "window of 25*bin"))

ggplot(ISdf_score, aes(x = V4, y = score, color = window)) +
  geom_line(linewidth = 1) +
  labs(x = "distance, kb", y = "insulation score") +
  scale_color_brewer(palette = "Dark2") +
  scale_x_continuous(guide = "prism_offset", limits = c(0, 60),
                     breaks = c(0, 20, 30, 40, 60),
                     labels = c("-15", "-5", "0", "5", "15")) +
  theme_prism(base_line_size = 0.5, base_family = "sans", base_fontface = "plain") +
  theme(text = element_text(size = 12),
        legend.title = element_blank(), legend.position = "top",
        legend.text = element_text(margin = margin(r = 10)),
        plot.margin = unit(c(0.8, 0.5, 0.5, 0), "cm")) +
  facet_grid(~name)
ggsave("mlei_IS_selection_EDF5a.pdf", width = 9, height = 3)


### Generate a final set of boundaries. For M. leidyi, it corresponds to boundaries
### called at the resolution 1000 bp with sliding window sizes 5 kb and 10 kb.
##########
df = read.csv("./mlei_1000bpRes.boundaries_Li.tsv", header = TRUE, sep = "\t")
df = subset(df, df$is_bad_bin=="False")

w1 = df[,c(1:3,7,13,16)]
colnames(w1) = c("chrom", "start", "end", "log2_IS", "boundary_strength", "is_boundary")
w1 = subset(w1, w1$is_boundary == "True")

w2 = df[,c(1:3,9,14,17)]
colnames(w2) = c("chrom", "start", "end", "log2_IS", "boundary_strength", "is_boundary")
w2 = subset(w2, w2$is_boundary=="True")

final = merge_w1_w2(w1, w2, "mlei_1000bp")