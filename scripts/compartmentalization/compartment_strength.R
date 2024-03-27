library("ggplot2")
library("reshape2")
library("plyr")
library("ggthemes")
library("ggprism")
library("dplyr")

################ 5K bins
################
setwd("./data/compartmentalization/5K_bins/")
filepath <- Sys.glob("*B-A.tsv")
c_score5K <- do.call(rbind, lapply(filepath, function(x) {
  df <- read.csv(x, sep = "\t", header = FALSE, stringsAsFactors = FALSE, row.names = NULL)
  df <- df[c(-1, -40), c(-1, -40)]
  BB <- sum(df[1:8, 1:8])
  AA <- sum(df[31:38, 31:38])
  AB <- sum(df[1:8, 31:38])
  BA <- sum(df[31:38, 1:8])
  total <- (AA + BB) / (AB + BA)
  totalA <- (AA) / (AB + BA)
  totalB <- (BB) / (AB + BA)
  total_df <- data.frame("AA" = totalA,
                         "BB" = totalB,
                         "species" = gsub("_.*", "", x),
                         "AA_percent" = totalA * 100 / total,
                         "BB_percent" = totalB * 100 / total)
  total_df
}))

c_score5K_m <- melt(c_score5K[,1:3])
names(c_score5K_m)[3] <- "5K"

label5K <- melt(c_score5K[,3:5])
label5K$variable <- gsub("_percent", "", label5K$variable)
names(label5K)[3] <- "5K" 
################

################ 10K bins
################
setwd("./data/compartmentalization/10K_bins/")
filepath <- Sys.glob("*B-A.tsv")
c_score10K <- do.call(rbind, lapply(filepath, function(x) {
  df <- read.csv(x, sep = "\t", header = FALSE, stringsAsFactors = FALSE, row.names = NULL)
  df <- df[c(-1, -40), c(-1, -40)]
  BB <- sum(df[1:8, 1:8])
  AA <- sum(df[31:38, 31:38])
  AB <- sum(df[1:8, 31:38])
  BA <- sum(df[31:38, 1:8])
  total <- (AA + BB) / (AB + BA)
  totalA <- (AA) / (AB + BA)
  totalB <- (BB) / (AB + BA)
  total_df <- data.frame("AA" = totalA,
                         "BB" = totalB,
                         "species" = gsub("_.*", "", x),
                         "AA_percent" = totalA * 100 / total,
                         "BB_percent" = totalB * 100 / total)
  total_df
}))

c_score10K_m <- melt(c_score10K[,1:3])
names(c_score10K_m)[3] <- "10K"

label10K <- melt(c_score10K[,3:5])
label10K$variable <- gsub("_percent", "", label10K$variable)
names(label10K)[3] <- "10K"
################

################ 20K bins
################
setwd("./data/compartmentalization/20K_bins/")
filepath <- Sys.glob("*B-A.tsv")
c_score20K <- do.call(rbind, lapply(filepath, function(x) {
  df <- read.csv(x, sep = "\t", header = FALSE, stringsAsFactors = FALSE, row.names = NULL)
  df <- df[c(-1, -40), c(-1, -40)]
  BB <- sum(df[1:8, 1:8])
  AA <- sum(df[31:38, 31:38])
  AB <- sum(df[1:8, 31:38])
  BA <- sum(df[31:38, 1:8])
  total <- (AA + BB) / (AB + BA)
  totalA <- (AA) / (AB + BA)
  totalB <- (BB) / (AB + BA)
  total_df <- data.frame("AA" = totalA,
                         "BB" = totalB,
                         "species" = gsub("_.*", "", x),
                         "AA_percent" = totalA * 100 / total,
                         "BB_percent" = totalB * 100 / total)
  total_df
}))

c_score20K_m <- melt(c_score20K[,1:3])
names(c_score20K_m)[3] <- "20K"

label20K <- melt(c_score20K[,3:5])
label20K$variable <- gsub("_percent", "", label20K$variable)
names(label20K)[3] <- "20K"
################


################ 50K bins
################
setwd("./data/compartmentalization/50K_bins/")
filepath <- Sys.glob("*B-A.tsv")
c_score50K <- do.call(rbind, lapply(filepath, function(x) {
  df <- read.csv(x, sep = "\t", header = FALSE, stringsAsFactors = FALSE, row.names = NULL)
  df <- df[c(-1, -40), c(-1, -40)]
  BB <- sum(df[1:8, 1:8])
  AA <- sum(df[31:38, 31:38])
  AB <- sum(df[1:8, 31:38])
  BA <- sum(df[31:38, 1:8])
  total <- (AA + BB) / (AB + BA)
  totalA <- (AA) / (AB + BA)
  totalB <- (BB) / (AB + BA)
  total_df <- data.frame("AA" = totalA,
                         "BB" = totalB,
                         "species" = gsub("_.*", "", x),
                         "AA_percent" = totalA * 100 / total,
                         "BB_percent" = totalB * 100 / total)
  total_df
}))

c_score50K_m <- melt(c_score50K[,1:3])
names(c_score50K_m)[3] <- "50K"

label50K <- melt(c_score50K[,3:5])
label50K$variable <- gsub("_percent", "", label50K$variable)
names(label50K)[3] <- "50K" 
################

## Join all dataframes
p1 <- join_all(list(c_score5K_m, c_score10K_m, c_score20K_m, c_score50K_m),
               by = c("species", "variable"), type = "full")

p1_sum <- aggregate(.~ species, p1[,-2], sum)
p1_sum_melt <- melt(p1_sum, id.vars = "species")

p1_sum_melt$species <- factor(p1_sum_melt$species,
                              levels = c("cowc", "mlei", "emue", "tadh", "dmel", "hsap"),
                              labels = c("Cowc", "Mlei", "Emue", "Tadh", "Dmel", "Hsap"))
p1_sum_melt$variable <- factor(p1_sum_melt$variable,
                               levels = c("5K", "10K", "20K", "50K"),
                               labels = c("5,000 bins", "10,000 bins", "20,000 bins", "50,000 bins"))

## Calculate standart deviation error bars
bars <- aggregate(. ~ species, p1_sum_melt[,-2], function(x) c(mean = mean(x), sd = sd(x)))
bars_value <- unlist(bars$value) 
error_bars <-data.frame("species" = bars$species, "mean" = bars_value[,1], "sd" = bars_value[,2])

## Plot compartment strength
ggplot(p1_sum_melt, aes(x=species, y=value, color = variable)) +
  geom_jitter(position=position_jitter(w=0.1, h=0.2), size = 2) +
  geom_point(data=error_bars, aes(x=species, y=mean), colour="black") +
  geom_errorbar(data = error_bars, 
                mapping = aes(x = species, y = mean, ymin = mean - sd, ymax = mean + sd),
                colour="black", width=.1) +
  scale_color_grey() +
  scale_x_discrete(guide = "prism_offset") +
  scale_y_continuous(guide = "prism_offset", limits = c(1.0, 2.5)) +
  labs(x = "", y = "Compartment strength\n(AA+BB) / (AB+BA)") +
  theme_prism(base_line_size = 0.5, base_family = "sans", base_fontface = "plain") +
  theme(axis.line.y = element_line(),
        text = element_text(size = 12), axis.text.x = element_text(angle = 30),
        legend.title = element_blank(), legend.text = element_text(margin = margin(r = 10)),
        plot.margin = unit(c(0.5,0.5,0,0), "cm")) 

ggsave("./compartment_strength_Fig2b.pdf", width = 5, height = 2.5)



## A and B ratio
################
ratio <- join_all(list(label5K, label10K, label20K, label50K),
                  by = c("species", "variable"), type = "full")
ratio <- melt(ratio, id.vars = c("species", "variable"))
ratio_ag <- aggregate(value ~ ., ratio[,-3], mean)
ratio_ag$species <- factor(ratio_ag$species,
                           levels = c("cowc", "mlei", "emue", "tadh", "dmel", "hsap"),
                           labels = c("Cowc", "Mlei", "Emue", "Tadh", "Dmel", "Hsap"))

ggplot(ratio_ag, aes(x=species, y=value,fill = variable)) +
  geom_bar(position="dodge", stat="identity") +
  geom_text(aes(label=round(value,0)), hjust=0.5, vjust=-1.4, 
                colour="black", size = 6, angle=0, position = position_dodge(.9)) +
  scale_fill_manual(values = c("#f56652", "#0060d2")) +
  labs(x = "", y = "% of total") +
  scale_y_continuous(guide = "prism_offset", trans = "reverse", breaks = c(0,60),labels = c(0,60)) +
  theme_prism(base_line_size = 0.5, base_family = "sans", base_fontface = "plain") +
  theme(axis.line.x = element_blank(), axis.line.y = element_line(), axis.ticks.x = element_blank(),
        text = element_text(size = 12), axis.text.x = element_blank(),
        legend.title = element_blank(), legend.text = element_text(margin = margin(r = 10)), legend.position = "bottom",
        plot.margin = unit(c(0.5,0.5,0,0), "cm")) +
  guides(fill = guide_legend(reverse = FALSE))

ggsave("./Fig2/AA_BB_percentage_Fig2b.pdf", width = 5, height = 2)