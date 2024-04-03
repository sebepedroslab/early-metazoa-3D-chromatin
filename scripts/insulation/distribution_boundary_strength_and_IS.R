library("rtracklayer")
library("GenomicRanges")
library("ggplot2")
library("ggforce")
library("ggridges")
library("ggprism")

### Functions
#############
w_and_s_boundary_insulation <- function(input_fn, input_li, input_li_is) {
  # insulation boundaries defined by applied Li threshold as strong
  names(mcols(input_li_is)) <- "insulation_score"
  names(mcols(input_li)) <- "boundary_strength"
  strong_ov <- findOverlaps(input_li, input_li_is)
  strong <- input_li[queryHits(strong_ov)]
  strong$insulation_score <- input_li_is$insulation_score[subjectHits(strong_ov)]

  # all annotated boundaries (weak and strong)
  input <- read.csv(input_fn, h = TRUE, sep = "\t")
  input <- input[,c(1:3,7,9,11:14)]
  colnames(input) <- c("seqnames", "start", "end",
                       "log2_insulation_score_w1", "log2_insulation_score_w2",
                       "boundary_strength_w1", "boundary_strength_w2",
                       "is_boundary_w1", "is_boundary_w2")
  input$coord <- paste(input$seqnames, input$start, input$end, sep = ":")
  input <- subset(input, input$is_boundary_w1 == "True" | input$is_boundary_w2 == "True")
  input$insulation_score <- ifelse(is.na(input$log2_insulation_score_w2), input$log2_insulation_score_w1, input$log2_insulation_score_w2)
  input$boundary_strength <- ifelse(is.na(input$boundary_strength_w2), input$boundary_strength_w1, input$boundary_strength_w2)
  input_gr <- makeGRangesFromDataFrame(input[,c(1:3,10,11,12)], keep.extra.columns = TRUE)

  # Define only weak boundaries
  ov <- findOverlaps(input_gr, strong)
  ovdf <- unique(data.frame("coord" = input_gr$coord[queryHits(ov)], strong[subjectHits(ov)]))
  weak_bound <- subset(input, !(input$coord %in% ovdf$coord))
  weak_bound <- weak_bound[,c(1:3,11,12)]

  # Label weak and strong insulation boundaries
  weak_bound$color <- "weak"
  strong_bound <- ovdf[,c(2:4,7:8)]
  strong_bound$color <- "strong"
  all_region <- rbind(weak_bound, strong_bound)
  all_region
}
#############

setwd("./data/insulation")

### Input insulation boundaries (weak and strong) from Step 4.
#######
cowc_fn <- "./cowc_400bpRes.boundaries_w+s.tsv"
emue_fn <- "./emue_1000bpRes.boundaries_w+s.tsv"
mlei_fn <- "./mlei_1000bpRes.boundaries_w+s.tsv"
tadh_fn <- "./tadh_400bpRes.boundaries_w+s.tsv"
dmel_fn <- "./dmel_400bpRes.boundaries_w+s.tsv"
hsap_fn <- "./hsap_10kbRes.boundaries_w+s.tsv"
#######

### Strong boundaries selected using Li threshold
######
cowc_li <- import("./boundary_strength/cowc_400bp_boundary_strength.bw")
cowc_li_is <- import("./insulation_score/cowc_400bp_insulation_score.bw")

emue_li <- import("./boundary_strength/emue_1000bp_boundary_strength.bw")
emue_li_is <- import("./insulation_score/emue_1000bp_insulation_score.bw")

mlei_li <- import("./boundary_strength/mlei_1000bp_boundary_strength.bw")
mlei_li_is <- import("./insulation_score/mlei_1000bp_insulation_score.bw")

tadh_li <- import("./boundary_strength/tadh_400bp_boundary_strength.bw")
tadh_li_is <- import("./insulation_score/tadh_400bp_insulation_score.bw")

dmel_li <- import("./boundary_strength/dmel_400bp_boundary_strenth.bw")
dmel_li_is <- import("./insulation_score/dmel_400bp_insulation_score.bw")

hsap_li <- import("./boundary_strength/hsap_10kb_boundary_strength.bw")
hsap_li_is <- import("./insulation_score/hsap_10kb_insulation_score.bw")
#######

### Identify weak and strong boundaries
#######
cowc <- w_and_s_boundary_insulation(cowc_fn, cowc_li, cowc_li_is)
emue <- w_and_s_boundary_insulation(emue_fn, emue_li, emue_li_is)
mlei <- w_and_s_boundary_insulation(mlei_fn, mlei_li, mlei_li_is)
tadh <- w_and_s_boundary_insulation(tadh_fn, tadh_li, tadh_li_is)
dmel <- w_and_s_boundary_insulation(dmel_fn, dmel_li, dmel_li_is)
hsap <- w_and_s_boundary_insulation(hsap_fn, hsap_li, hsap_li_is)

all <- data.frame("species" = c(rep("Cowc", length(cowc$boundary_strength)),
                                rep("Emue", length(emue$boundary_strength)),
                                rep("Mlei", length(mlei$boundary_strength)),
                                rep("Tadh", length(tadh$boundary_strength)),
                                rep("Dmel", length(dmel$boundary_strength)),
                                rep("Hsap", length(hsap$boundary_strength))),
                  "insulation" = c(cowc$insulation_score,
                                   emue$insulation_score,
                                   mlei$insulation_score,
                                   tadh$insulation_score,
                                   dmel$insulation_score,
                                   hsap$insulation_score),
                  "boundary" = c(cowc$boundary_strength,
                                 emue$boundary_strength,
                                 mlei$boundary_strength,
                                 tadh$boundary_strength,
                                 dmel$boundary_strength,
                                 hsap$boundary_strength),
                  "category" = c(cowc$color,
                                 emue$color,
                                 mlei$color,
                                 tadh$color,
                                 dmel$color,
                                 hsap$color))

all$species <- factor(all$species, levels = c("Cowc", "Mlei", "Emue","Tadh","Dmel", "Hsap"))

# Prepare a dataframe with a umber of final strong and weak boundaries
all_label <- all[,-2]
all_label$boundary <- 1
all_label <- aggregate(boundary ~ ., all_label, sum)
all_label$total_number <- paste(all_label$category, " = ", all_label$boundary, sep = "")
all_label <- all_label[,c(1,4)]
all_label <- aggregate(total_number ~ species, all_label, paste, collapse = "\n")
all_label$weak <- all_label$total_number
all_label$weak <- gsub(".*\n", "", all_label$weak)
all_label$strong <- gsub("\n.*", "", all_label$total_number)

# Plot the distribution of insulation boundary strength values
ggplot(all, aes(x = species, y = boundary)) +
  geom_violin() +
  geom_sina(aes(colour = factor(category), group = factor(species)), size = 0.1, alpha = 0.4, bins = 1) +
  geom_boxplot(width = 0.05, outlier.shape = NA, linewidth = 0.2) +
  labs(x = "", y = "boundary strength") +
  scale_color_manual(values = c("#5c8ca3ff", "grey")) +
  geom_text(data = all_label ,aes(x = c(0.8, 1.8, 2.8, 3.8, 4.8, 5.8), y = 4, label = strong),
            color = "#5c8ca3ff", angle = 90, size = 3, hjust = 0) +
  geom_text(data = all_label ,aes(x = c(1.2, 2.2, 3.2, 4.2, 5.2, 6.2), y = 4, label = weak),
            color = "grey30", angle = 90, size = 3, hjust = 0) +
  scale_y_continuous(trans = "log10", guide = "prism_offset", limits = c(0.00001, 100),
                     labels = function(x) format(round(x,3), scientific = TRUE),
                     breaks = c(0.00001, 0.01, 10)) +
  theme_prism(base_line_size = 0.5, base_family = "sans", base_fontface = "plain") +
  theme(axis.line.y = element_line(), axis.line.x = element_blank(),
        text = element_text(size = 12), legend.title = element_blank(),
        axis.text.x = element_text(angle = 0), axis.ticks.x = element_blank(),
        legend.text = element_text(margin = margin(r = 10)), legend.position = "none",
        plot.margin = unit(c(0.8,0.5,0.5,0), "cm"))
ggsave("./All_species_boundary_strength_EDF5c.pdf", width = 4.8, height = 4.8)

# Plot the distribution of the insulation score values
all$category <- gsub("strong", "strong boundaries", all$category)
all$category <- gsub("weak", "weak boundaries", all$category)

ggplot(all, aes(x = insulation, y = species, fill = category)) +
  geom_density_ridges2(scale = 1, quantile_lines = TRUE, quantiles = 0.5, size = 0.3) +
  scale_x_continuous(guide = "prism_offset", limits = c(-2,0.5)) +
  labs(x = "Insulation score", y = "") +
  scale_fill_manual(values = c("#5c8ca3ff", "#bebebe33")) +
  theme_prism(base_line_size = 0.5, base_family = "sans", base_fontface = "plain") +
  theme(text = element_text(size = 12),
        legend.title = element_blank(), legend.position = "top",
        legend.text = element_text(margin = margin(r = 10)),
        plot.margin = unit(c(0.8,0.5,0.5,0), "cm"),
        axis.line.y = element_blank(), axis.ticks.y = element_blank())
ggsave("./All_species_density_insulation_Fig2d.pdf", width = 4, height = 3.5)