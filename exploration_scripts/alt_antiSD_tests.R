
# This is a script used for testing is a "mixed" 16S rRNA gene set is uniformly distributed accross
#   the bacterial phyla.

phyla <- c('Proteobacteria', 'Firmicutes', 'Actinobacteria', 'Bacteroidetes',
           'Tenericutes', 'Cyanobacteria', 'Chlamydiae', 'Spirochaetes','Verrucomicrobia')
num_genomes <- c(14311, 5372, 2444, 905, 423, 188, 181, 176, 112)


# "Mixed" 16S rRNA gene sets

antiSD_counts <- c(31, 19, 10, 1, 1, 0, 0, 0, 0)
mid_no_antiSD_df <- data.frame(
  phylum = phyla,
  num_genomes = num_genomes,
  antiSD_counts = antiSD_counts,
  noantiSD_counts = num_genomes - antiSD_counts
)

fisher.test(as.matrix(mid_no_antiSD_df[, 3:4]), simulate.p.value=F)

