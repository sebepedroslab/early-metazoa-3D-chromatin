library(ggplot2)

setwd("./data/chromatin_loops/SIP")
score_quantiles <- c(0.00, 0.025, 0.05, 0.10, 0.15, 0.25, 0.50, 0.75, 0.80, 0.85, 0.90, 0.95, 0.97, 0.98, 0.99, 0.995, 1.00)

### For Mnemiopsis leidyi
#########################
loops <- read.csv("./Mlei_finalLoops.txt", header = TRUE, sep = "\t")
value_quantiles <- quantile(loops$APScoreAvg, score_quantiles)

# Define the cutoff value based on the peak prominence of APScoreAvg values
ggplot(loops, aes(x= APScoreAvg)) +
  geom_density() +
  scale_x_continuous(trans = "log10") +
  geom_vline(xintercept = value_quantiles["99.5%"], colour = "red", linetype = "dashed")

aps_loops <- subset(loops, loops$APScoreAvg < value_quantiles["99.5%"]) 
write.table(aps_loops, "./data/chromatin_loops/Mlei_filtered_chromatin_loops_4261.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
#########################


### For Trichoplax adhaerens
############################
loops <- read.csv("./Tadh_finalLoops.txt", header = TRUE, sep = "\t")
value_quantiles <- quantile(loops$APScoreAvg, score_quantiles)

ggplot(loops, aes(x= APScoreAvg)) +
  geom_density() +
  scale_x_continuous(trans = "log10") +
  geom_vline(xintercept = value_quantiles["97%"], colour = "red", linetype = "dashed")

aps_loops = subset(loops, loops$APScoreAvg < value_quantiles["97%"])
write.table(aps_loops, "./data/chromatin_loops/Tadh_filtered_chromatin_loops_3065.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
############################