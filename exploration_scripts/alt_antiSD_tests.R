
# This is a script used for testing is a "mixed" 16S rRNA gene set is uniformly distributed accross
#   the bacterial phyla.

phyla <- c('Proteobacteria', 'Firmicutes', 'Actinobacteria', 'Bacteroidetes',
           'Tenericutes', 'Cyanobacteria', 'Chlamydiae', 'Spirochaetes','Verrucomicrobia')
num_genomes <- c(13144, 4976, 2260, 816, 406, 184, 174, 157, 110)


# "Mixed" 16S rRNA gene sets

antiSD_counts <- c(30, 18, 9, 1, 1, 0, 0, 0, 0)
mid_no_antiSD_df <- data.frame(
  phylum = phyla,
  num_genomes = num_genomes,
  antiSD_counts = antiSD_counts,
  noantiSD_counts = num_genomes - antiSD_counts
)

fisher.test(as.matrix(mid_no_antiSD_df[, 3:4]), simulate.p.value=TRUE)
