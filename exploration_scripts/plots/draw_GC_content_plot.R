library(ggplot2)


infpath <- '/mnt/1.5_drive_0/RiboGrove_workdirs/4.210/per_genome_avg_16S_gc_df.tsv'
df <- read.csv(infpath, sep = '\t')

# df <- subset(df, df$domain == 'Bacteria')
df$gc_16S <- df$gc_16S * 100.0
df$gc_WG <- df$gc_WG * 100.0

# abund_phyla <- c('Proteobacteria', 'Firmicutes', 'Actinobacteria', 'Bacteroidetes')
abund_phyla <- c('Proteobacteria', 'Firmicutes', 'Actinobacteria')
df$Organism <- df$phylum

for (i in 1:nrow(df)) {
  if (! (df$phylum[i] %in% abund_phyla)) {
    df$Organism[i] <- 'Other bacteria\nand archaea'
  }
}

proteobacteria_color <- '#000000'#'#8B58C1'
firmicutes_color <- '#5A5A5A'
actinobacteria_color <- '#A3A3A3'
# bacteroidetes_color <- '#6CBB3C'
other_color <- '#BFBFBF'

proteobacteria_shape <- 20
firmicutes_shape <- 16
actinobacteria_shape <- 17
# bacteroidetes_shape <- 18
other_shape <- 3

proteobacteria_alpha <- 1.0
firmicutes_alpha <- 0.4
actinobacteria_alpha <- 0.4
# bacteroidetes_alpha <- 18
other_alpha <- 0.6

proteobacteria_size <- 0.4
firmicutes_size <- 1.7
actinobacteria_size <- 1.9
# bacteroidetes_size <- 18
other_size <- 1.4

# colors <- c(actinobacteria_color, bacteroidetes_color, firmicutes_color, other_color, proteobacteria_color)
# shapes <- c(actinobacteria_shape, bacteroidetes_shape, firmicutes_shape, other_shape, proteobacteria_shape)
# alphas <- c(0.4, 0.4, 0.4, 0.4, 0.5)
# sizes <- c(1.9, 2.3, 1.7, 1.5, 1.5)
colors <- c(actinobacteria_color, firmicutes_color, other_color, proteobacteria_color)
shapes <- c(actinobacteria_shape, firmicutes_shape, other_shape, proteobacteria_shape)
alphas <- c(actinobacteria_alpha, firmicutes_alpha, other_alpha, proteobacteria_alpha)
sizes <- c(actinobacteria_size, firmicutes_size, other_size, proteobacteria_size)


ggplot(data = subset(df, df$Organism != 'Proteobacteria'),
       aes(x=gc_WG, y = gc_16S, col = Organism, alpha = Organism, shape = Organism, size = Organism)) + 
  geom_point() +
  geom_point(data = subset(df, df$Organism == 'Proteobacteria')) +
  scale_color_manual(values=colors) +
  scale_alpha_manual(values=alphas) +
  scale_shape_manual(values=shapes) +
  scale_size_manual(values=sizes) +
  labs(x = 'Whole genome GC content, (%)', y = 'Average 16S rRNA GC content, (%)') +
  theme_bw()



# == Color ==

abund_phyla <- c('Proteobacteria', 'Firmicutes', 'Actinobacteria', 'Bacteroidetes')
# abund_phyla <- c('Proteobacteria', 'Firmicutes', 'Actinobacteria')
df$Organism <- df$phylum

for (i in 1:nrow(df)) {
  if (! (df$phylum[i] %in% abund_phyla)) {
    df$Organism[i] <- 'Other bacteria\nand archaea'
  }
}

proteobacteria_color <- '#000000'#'#8B58C1'
firmicutes_color <- '#4AA0F0'
actinobacteria_color <- '#FFD257'
bacteroidetes_color <- '#5EBA4C'
other_color <- '#9F9F9F'

proteobacteria_shape <- 20
firmicutes_shape <- 15
actinobacteria_shape <- 17
bacteroidetes_shape <- 18
other_shape <- 4

proteobacteria_alpha <- 1.0
firmicutes_alpha <- 0.5
actinobacteria_alpha <- 0.5
bacteroidetes_alpha <- 0.5
other_alpha <- 0.5

proteobacteria_size <- 0.4
firmicutes_size <- 1.7
actinobacteria_size <- 1.9
bacteroidetes_size <- 2.3
other_size <- 1.0

# colors <- c(actinobacteria_color, bacteroidetes_color, firmicutes_color, other_color, proteobacteria_color)
# shapes <- c(actinobacteria_shape, bacteroidetes_shape, firmicutes_shape, other_shape, proteobacteria_shape)
# alphas <- c(0.4, 0.4, 0.4, 0.4, 0.5)
# sizes <- c(1.9, 2.3, 1.7, 1.5, 1.5)
colors <- c(actinobacteria_color, bacteroidetes_color, firmicutes_color, other_color, proteobacteria_color)
shapes <- c(actinobacteria_shape, bacteroidetes_shape, firmicutes_shape, other_shape, proteobacteria_shape)
alphas <- c(actinobacteria_alpha, bacteroidetes_alpha, firmicutes_alpha, other_alpha, proteobacteria_alpha)
sizes <- c(actinobacteria_size, bacteroidetes_size, firmicutes_size, other_size, proteobacteria_size)


ggplot(data = subset(df, df$Organism != 'Proteobacteria'),
       aes(x=gc_WG, y = gc_16S, col = Organism, alpha = Organism, shape = Organism, size = Organism)) + 
  geom_point() +
  geom_point(data = subset(df, df$Organism == 'Proteobacteria')) +
  scale_color_manual(values=colors) +
  scale_alpha_manual(values=alphas) +
  scale_shape_manual(values=shapes) +
  scale_size_manual(values=sizes) +
  labs(x = 'Whole genome GC content, (%)', y = 'Average 16S rRNA GC content, (%)') +
  theme_bw()
