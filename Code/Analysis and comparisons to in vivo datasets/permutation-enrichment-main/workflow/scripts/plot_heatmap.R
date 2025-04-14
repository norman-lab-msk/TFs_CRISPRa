# This script plots a heatmap of the permutation testing results.
#
# Author: Karthik Guruvayurappan

library(dplyr)
library(ggplot2)

# read in the dataframes containing permutation testing results
results <- data.frame()

for (i in 1:32) {
    batch_results <- read.csv(snakemake@input[[i]])
    results <- rbind(results, batch_results)
}

# for (i in 1:32) {
#   batch_results <- read.csv(paste0('results/permutation_results/permutation_results_', i, '.csv'))
#   results <- rbind(results, batch_results)
# }
# modify data types
results$Program <- as.factor(results$Program)
results$adj_p <- p.adjust(results$pval, method = 'fdr')

# write to output CSV
write.csv(
  results,
  snakemake@output[[1]],
  row.names = FALSE,
  quote = FALSE
)

# write.csv(
#   results,
#   'results/permutation_results.csv',
#   row.names = FALSE,
#   quote = FALSE
# )

# add binary significant column
results$significant <- as.integer(results$adj_p < 0.1)

# programs to plot (only plot programs with overlaps)
program_counts <- results %>% group_by(Program) %>% summarize(total = sum(significant))
plot_programs <- (program_counts[program_counts$total >= 1, ]$Program)
results <- results[results$Program %in% plot_programs, ]

# change significance to factor
results$significant <- as.factor(results$significant)

plot <- ggplot(results, aes(Program, Cluster, fill = significant, alpha = overlap_coefs)) +
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
  filename = snakemake@output[[2]],
  unit = 'in',
  height = 12,
  width = 9
)
