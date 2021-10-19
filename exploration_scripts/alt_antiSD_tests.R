phyla <- c('Proteobacteria', 'Firmicutes', 'Actinobacteria', 'Bacteroidetes',
          'Tenericutes', 'Spirochaetes', 'Chlamydiae', 'Cyanobacteria', 'Verrucomicrobia')
# num_genomes <- c(14859, 5576, 2520, 948, 427, 256, 244, 233, 116)

phyla <- c('Proteobacteria', 'Firmicutes', 'Actinobacteria', 'Bacteroidetes',
           'Tenericutes', 'Cyanobacteria', 'Chlamydiae', 'Spirochaetes','Verrucomicrobia')
num_genomes <- c(13144, 4976, 2260, 816, 406, 184, 174, 157, 110)


# Repeats

num_repeat_genes <- c(32, 14, 1, 0, 0, 0, 0, 0, 0)

rep_df <- data.frame(
  phylum = phyla,
  num_genomes = num_genomes,
  num_repeat_genes = num_repeat_genes,
  num_nonrep_genes = num_genomes - num_repeat_genes
)

# chisq.test(as.matrix(rep_df[, 3:4]))
# fisher.test(as.matrix(rep_df[, 3:4]))

for (i in 3:length(phyla)) {
  p <- fisher.test(as.matrix(rep_df[1:i, 3:4]))$p.value
  cat(sprintf("%d: %f\n", i, p))
}


# Low-GC IVSs

# num_lowGCivs_genes <- c(4, 12, 0, 0, 0, 0, 0, 0, 0)
num_lowGCivs_genes <- c(4, 8, 0, 0, 0, 0, 0, 0, 0)
lowGCivs_df <- data.frame(
  phylum = phyla,
  num_genomes = num_genomes,
  num_lowGCivs_genes = num_lowGCivs_genes,
  num_normal_genes = num_genomes - num_lowGCivs_genes
)

# 
# chisq.test(as.matrix(lowGCivs_df[, 3:4]))
# fisher.test(as.matrix(lowGCivs_df[, 3:4]))

for (i in 2:length(phyla)) {
  p <- fisher.test(as.matrix(lowGCivs_df[1:i, 3:4]))$p.value
  cat(sprintf("%d: %f\n", i, p))
}



# Zero-antiSD DF

# phyla <- c('Actinobacteria', 'Bacteroidetes', 'Chlamydiae', 'Chloroflexi',
#            'Firmicutes', 'Proteobacteria', 'Tenericutes')
# num_genomes <- c(2520, 948, 244, 45, 5576, 14859, 427)
# no_antiSD_counts <- c(1, 244, 1, 2, 3, 25, 8)

no_antiSD_counts <- c(25, 3, 1, 244, 8, 0, 1, 0, 0)
zero_no_antiSD_df <- data.frame(
  phylum = phyla,
  num_genomes = num_genomes,
  no_antiSD_counts = no_antiSD_counts,
  antiSD_counts = num_genomes - no_antiSD_counts
)

fisher.test(as.matrix(zero_no_antiSD_df[, 3:4]), simulate.p.value=TRUE)


# Mid-antiSD DF

# phyla <- c('Actinobacteria', 'Bacteroidetes', 'Firmicutes',
#            'Proteobacteria', 'Synergistetes', 'Tenericutes')
# num_genomes <- c(2520, 948, 5576, 14859, 5, 427)
# antiSD_counts <- c(9, 1, 23, 36, 1, 1)


antiSD_counts <- c(30, 18, 9, 1, 1, 0, 0, 0, 0)
mid_no_antiSD_df <- data.frame(
  phylum = phyla,
  num_genomes = num_genomes,
  antiSD_counts = antiSD_counts,
  noantiSD_counts = num_genomes - antiSD_counts
)

fisher.test(as.matrix(mid_no_antiSD_df[, 3:4]), simulate.p.value=TRUE)

sum(mid_no_antiSD_df$antiSD_counts)


