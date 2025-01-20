library("qdap")
library("ggplot2")
library("mclust")
library("mixtools")
library("scales")

### Functions
##############
plot_mix_comps <- function(x, mu, sigma, lam) {
  lam * dnorm(x, mu, sigma)
}

eigen_assignment <- function(input_fn, species, width_chr, height_chr) {
  df_e123 <- read.csv(input_fn, h = TRUE, sep = "\t")
  df_e123 <- df_e123[,c(1:4)]
  df_e123 <- na.omit(df_e123)
  q1 <- quantile(df_e123$E1, 0.025)
  q2 <- quantile(df_e123$E1, 0.975)
  df <- subset(df_e123, df_e123$E1 > q1 & df_e123$E1 < q2)
  df$E1_norm <- df$E1

  pdf(paste(species, "eigenvalues_20K_clusters.pdf", sep = "_"), width = 3, height = 3.2)
  plot(mclustBIC(df$E1_norm))  ### identify number of clusters
  dev.off()

  w1 <- normalmixEM(df$E1_norm, k = 3)
  plot_df <- data.frame("V1" = w1$x)

  z <- seq(min(df$E1_norm), max(df$E1_norm), 1e-2)
  u <- sapply(1:length(w1$mu), FUN = function(i) dnorm(z, w1$mu[i], w1$sigma[i]))

  # find intersection points by retrieving the values with minimal absolute difference of the likelihoods given neighboring distributions
  inflections <- sapply(1:(length(w1$mu) - 1), FUN = function(i){
    z_range <- z[z > w1$mu[i] & z < w1$mu[i + 1]]
    u_range <- u[z > w1$mu[i] & z < w1$mu[i + 1],]
    inf <- z_range[which.min(abs(u_range[,i] * w1$lambda[i] - u_range[,i + 1] * w1$lambda[i + 1]))]
    return(inf)
  })

  ggplot(plot_df, aes(V1)) +
    geom_density(fill = "#d9d9d9", trim = TRUE) +
    stat_function(geom = "line", fun = plot_mix_comps,
                  args = list(w1$mu[1], w1$sigma[1], lam = w1$lambda[1]),
                  colour = "#0060d2") +
    stat_function(geom = "line", fun = plot_mix_comps,
                  args = list(w1$mu[2], w1$sigma[2], lam = w1$lambda[2]),
                  colour = "#fed976") +
    stat_function(geom = "line", fun = plot_mix_comps,
                  args = list(w1$mu[3], w1$sigma[3], lam = w1$lambda[3]),
                  colour = "#f56652") +
    labs(x = "Eigenvector E1", y = "Density", title = species) +
    scale_x_continuous(limits = c(-1, 1)) +
    scale_y_continuous(limits = c(0, 1.4)) +
    geom_vline(data = data.frame(x = inflections), aes(xintercept = x), col = "red", linetype = "dashed") +
    theme_bw(base_line_size = 0.5, base_family = "sans") +
    theme(text = element_text(size = 12), panel.grid = element_blank(), plot.title = element_text(hjust = 0.5))
  ggsave(paste(species, "eigenvalues_gaus_20K_EDF4c.pdf", sep = "_"), width = 2.5, height = 1.8)

  df$comp <- ifelse(df$E1_norm <= inflections[1], "B",
                    ifelse(df$E1_norm >= inflections[2], "A", "I"))

  write.table(df[,c(1:3,6)], paste(species, "compartment_regions_Gaus_distr.txt", sep = "_"),
              col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

  ggplot(df, aes(x = E1_norm)) +
    geom_density(trim = TRUE) +
    labs(x = "Eigenvector E1", y = "Density", title = species) +
    scale_x_continuous(limits = c(-1, 1)) +
    scale_y_continuous(limits = c(0, 2)) +
    facet_wrap(~chrom, ncol = 11) +
    geom_vline(data = data.frame(x = inflections), aes(xintercept = x), col = "red", linetype = "dashed") +
    theme_bw(base_line_size = 0.5, base_family = "sans") +
    theme(text = element_text(size = 12), panel.grid = element_blank(), plot.title = element_text(hjust = 0.5))
  ggsave(paste(species, "eigen_chr_grid_Gaus_distribution_EDF4d.pdf", sep = "_"), width = width_chr, height = height_chr)

  df$species <- species
  df$count <- 1
  df
}
##############

setwd("./data/compartmentalization/eigenvalues")

## Path to the files
dmel_fn <- "./dmel_E1_6878bp.bed"
hsap_fn <- "./hsap_E1_154413bp.bed"
tadh_fn <- "./tadh_E1_4509bp.bed"
mlei_fn <- "./mlei_E1_10252bp.bed"
emue_fn <- "./emue_E1_10790bp.bed"
cowc_fn <- "./cowc_E1_1423bp.bed"
sros_fn <- "./sros_E1_2585bp.bed"
sarc_fn <- "./sarc_E1_7032bp.bed"
nvec_fn <- "./nvec_E1_13469bp.bed"

## Assign compartment types
dmel <- eigen_assignment(dmel_fn, "Dmel", 10, 3.5)
hsap <- eigen_assignment(hsap_fn, "Hsap", 10, 10)
tadh <- eigen_assignment(tadh_fn, "Tadh", 10, 3.5)
mlei <- eigen_assignment(mlei_fn, "Mlei", 10, 6.5)
emue <- eigen_assignment(emue_fn, "Emue", 10, 10)
cowc <- eigen_assignment(cowc_fn, "Cowc", 10, 10)
sros <- eigen_assignment(sros_fn, "Sros", 10, 10)
sarc <- eigen_assignment(sarc_fn, "Sarc", 10, 10)
nvec <- eigen_assignment(nvec_fn, "Nvec", 10, 3.5)

## Prepare data for stacked barplot
my_list <- list("sarc" = sarc, "cowc" = cowc, "sros" = sros, "emue" = emue,
                "mlei" = mlei, "tadh" = tadh, "nvec" = nvec, "dmel" = dmel, "hsap" = hsap)

do.call(rbind, lapply(my_list, function(x) {
  df <- x[,c(6,7,8)]
  df_agg <- aggregate(count ~ ., df, sum)
  df_agg$percent <- df_agg$count * 100 / sum(df_agg$count)
  df_agg$label <- paste(df_agg$comp, "\n", round(df_agg$percent), "%", sep = "")
  df_agg$comp <- factor(df_agg$comp, levels = c("A", "I", "B"))

  ggplot(df_agg, aes(x = species, y = percent, fill = comp, label = label)) +
    geom_bar(stat = "identity", colour = "black", linewidth = 0.05) +
    scale_fill_manual(values = c("#f56652", "#fed976",  "#0060d2")) +
    coord_flip() +
    geom_text(position = position_stack(vjust = 0.5)) +
    guides(fill = guide_legend(reverse = TRUE)) +
    theme_minimal(base_family = "sans") +
    theme(panel.grid = element_blank(), axis.text = element_blank(), axis.line = element_blank(), axis.title = element_blank(),
          text = element_text(size = 12), legend.title = element_blank(),
          legend.text = element_text(margin = margin(r = 10)), legend.position = "none",
          plot.margin = unit(c(0,0,0,0), "cm"))
  ggsave(paste(unique(df_agg$species), "Gaus_distribution_barplot_EDF3c.pdf", sep = "_"), width = 4, height = 0.7)
}))