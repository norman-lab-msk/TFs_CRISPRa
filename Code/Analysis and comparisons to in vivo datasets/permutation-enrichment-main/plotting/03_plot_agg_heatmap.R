# This script plots a heatmap of permutation testing results across all
# atlases.
#
# Author: Karthik Guruvayurappan

library(dplyr)
library(ggplot2)

# import results across all atlases
results_dir <- 'results/results_archive/'
atlases <- c(
    'buechler_atlas',
    'cardiac_fibro',
    'myofibroblast',
    'korsunsky_atlas'
)

results_df = data.frame()

for (i in 1:length(atlases)) {
    atlas_results <- read.csv(paste0(
        results_dir, atlases[i], '/permutation_results.csv'
    ))
    atlas_results$atlas <- atlases[i]
    results_df <- rbind(results_df, atlas_results)
}

# filter for programs of interest
keep_programs <- c(
    18, # universal
    31, # universal
    54, # ?
    27, # inflammatory
    32, # inflammatory
    1, # myofibroblast?
    4, # myofibroblast
    19 # immune
)

results_df <- results_df[results_df$Program %in% keep_programs, ]

# plot results across all clusters
plot_df <- results_df
plot_df$Program <- as.factor(plot_df$Program)
plot_df$significant <- as.integer(plot_df$adj_p < 0.1) # add binary significance

# clusters to plot (only plot clusters with overlaps)
cluster_counts <- plot_df %>% group_by(Cluster) %>% summarize(total = sum(significant))
plot_clusters <- (cluster_counts[cluster_counts$total >= 1, ]$Cluster)
plot_df <- plot_df[plot_df$Cluster %in% plot_clusters, ]

plot_df$significant <- as.factor(plot_df$significant)

plot <- ggplot(plot_df, aes(Program, Cluster, fill = significant, alpha = overlap_coefs)) +
  geom_tile(color = 'black') +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(color = 'black', size = 14),
    axis.text.y = element_text(color = 'black', size = 14),
    axis.title.x = element_text(color = 'black', size = 20),
    axis.title.y = element_text(color = 'black', size = 20)
  ) +
  scale_fill_manual(
    breaks = c('0', '1'),
    values = c('white', '#808080')
  ) +
  guides(fill = 'none') +
  scale_alpha_continuous(name = 'Overlap Coef', range = c(0, 1))

ggsave(
  plot = plot,
  device = 'pdf',
  filename = 'results/agg_heatmap_allclusters.pdf',
  unit = 'in',
  height = 12,
  width = 9
)

# keep specific clusters
keep_clusters <- c(
    # 'CD34+MFAP5+ C9 Synovium', # universal
    # 'CD34+MFAP5+ C9 SalivaryGland', # universal
    'CD34+MFAP5+ C9 Marginal', # universal
    # 'CD34+MFAP5+ C9 Lung', # universal
    # 'CD34+MFAP5+ C9 Gut', # universal
    'PI16+ (n = 779)', # universal
    'POLCE2.MFAP5 Fib', # universal
    'ground state Fib', # universal
    'c05', # universal
    'c03', # universal
    # 'SPARC+COL3A1+ C4 Synovium', # inflammatory
    # 'SPARC+COL3A1+ C4 SalivaryGland', # inflammatory
    'SPARC+COL3A1+ C4 Marginal', # inflammatory
    # 'SPARC+COL3A1+ C4 Lung', # inflammatory
    # 'SPARC+COL3A1+ C4 Gut', # inflammatory
    'COL3A1 (n = 72)', # inflammatory
    'LRRC15+ (n = 472)', # inflammatory
    'POSTN+ Fib', # inflammatory
    'c04', # inflammatory
    # 'MYH11+ C13 Synovium', # myofibroblast
    # 'MYH11+ C13 SalivaryGland', # myofibroblast
    'MYH11+ C13 Marginal', # myofibroblast
    # 'MYH11+ C13 Lung', # myofibroblast
    # 'MYH11+ C13 Gut', # myofibroblast
    'Myofibroblast', # myofibroblast
    'c01', # myofibroblast
    # 'CXCL10+CCL19+ C11 Synovium', # immune/antigen presentation 
    # 'CXCL10+CCL19+ C11 SalivaryGland', # immune
    'CXCL10+CCL19+ C11 Marginal', # immune
    # 'CXCL10+CCL19+ C11 Lung', # immune
    # 'CXCL10+CCL19+ C11 Gut', # immune
    'CCL19+ (n = 400)', # immune
    'Type I IFN Fib', # immune
    'c20' # immune
)

results_df <- results_df[results_df$Cluster %in% keep_clusters, ]
results_df$Cluster <- factor(results_df$Cluster, levels = keep_clusters)

# plot results across all clusters
plot_df <- results_df
plot_df$Program <- as.factor(plot_df$Program)
plot_df$significant <- as.integer(plot_df$adj_p < 0.1) # add binary significance

# clusters to plot (only plot clusters with overlaps)
cluster_counts <- plot_df %>% group_by(Cluster) %>% summarize(total = sum(significant))
plot_clusters <- (cluster_counts[cluster_counts$total >= 1, ]$Cluster)
plot_df <- plot_df[plot_df$Cluster %in% plot_clusters, ]

plot_df$significant <- as.factor(plot_df$significant)

plot <- ggplot(plot_df, aes(Program, Cluster, fill = significant, alpha = overlap_coefs)) +
  geom_tile(color = 'black') +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(color = 'black', size = 14),
    axis.text.y = element_text(color = 'black', size = 14),
    axis.title.x = element_text(color = 'black', size = 20),
    axis.title.y = element_text(color = 'black', size = 20)
  ) +
  scale_fill_manual(
    breaks = c('0', '1'),
    values = c('white', '#808080')
  ) +
  guides(fill = 'none') +
  scale_alpha_continuous(name = 'Overlap Coef', range = c(0, 1))

ggsave(
  plot = plot,
  device = 'pdf',
  filename = 'results/agg_heatmap_select_clusters.pdf',
  unit = 'in',
  height = 12,
  width = 9
)