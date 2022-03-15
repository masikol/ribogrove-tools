library(ggplot2)

infpath <- '/mnt/1.5_drive_0/RiboGrove_releases/4.210/metadata/gene_seqs_statistics.tsv'
stats_df <- read.csv(infpath, sep='\t')

stats_df['gc_skew'] <- (stats_df['g'] - stats_df['c']) /
                       (stats_df['g'] + stats_df['c'])

stats_df['at_skew'] <- (stats_df['a'] - stats_df['t']) /
                       (stats_df['a'] + stats_df['t'])

stats_df <- stats_df[complete.cases(stats_df[c('domain', 'phylum')]),]


set_Organism <- function(x, output) {
  result <- 'Other bacteria'
  if (x[['domain']] == 'Bacteria') {
    if (x[['phylum']] == 'Actinobacteria') {
      result <- 'Actinobacteria'
    }
  } else {
    result <- 'Archaea'
  }
  return(result)
}

stats_df$Organism <- as.character(apply(stats_df, 1, set_Organism))

# actinobacteria_color <- '#FFD445'
# actinobacteria_line_color <- '#FFC500'
# archaea_color <- '#50a1f5'
# archaea_line_color <- '#2A8BEF'
# otherbacteria_color <- '#7F7F7F'

# actinobacteria_shape <- 4
# archaea_shape <- 2
# otherbacteria_shape <- 19

# == Color ==

actinobacteria_color <- '#ffd258'
actinobacteria_line_color <- '#ffc337'
archaea_color <- '#48a0f1'
archaea_line_color <- '#208bea'
otherbacteria_color <- '#7e7e7e'

actinobacteria_shape <- 17
archaea_shape <- 15
otherbacteria_shape <- 21



ggplot(data = stats_df) +
  geom_point(
    data = subset(stats_df, stats_df$Organism == 'Other bacteria'),
    aes(
      x = gc_skew,
      y = at_skew,
      col = Organism,
      alpha = Organism,
      shape = Organism
    ),
    size=1.2
) +
  geom_point(
    data = subset(stats_df, stats_df$Organism != 'Other bacteria'),
    aes(
      x = gc_skew,
      y = at_skew,
      col = Organism,
      alpha = Organism,
      shape = Organism
    ),
    size=1.6
  ) +
  scale_color_manual(values=c(actinobacteria_color, archaea_color, otherbacteria_color)) +
  scale_alpha_manual(values=c(1.0, 1.0, 0.3)) +
  scale_shape_manual(values=c(actinobacteria_shape, archaea_shape, otherbacteria_shape)) +
  labs(x = 'GC skew', y = 'AT skew') +
  geom_smooth(
    data = subset(stats_df, stats_df$Organism == 'Actinobacteria'),
    aes(x = gc_skew, y = at_skew),
    method=lm,
    color = actinobacteria_line_color,
    se = F
    # linetype = actinobacteria_line_type
  )+
  geom_smooth(
    data = subset(stats_df, stats_df$Organism == 'Archaea'),
    aes(x = gc_skew, y = at_skew),
    method=lm,
    color = archaea_line_color,
    se = F
    # linetype = archaea_line_type
  ) +
  theme_bw()







