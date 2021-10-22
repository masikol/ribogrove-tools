Here, scripts used for exploration of RybaSom data are located.

# Content

1. `calculate_entropy.py` -- the script for calculating intragenomic veriability of 16S rRNA gene sequences. More specifically, it calculates Shannon entropy for each MSA position of sequences aligned.

2. `check_primers_mfeprimer.py` -- the script takes a fasta file of gene sequences and checks if hardcoded primers can produce some product with input sequences as templates.

3. `count_bases.py` -- the script counts all IUPAC bases in sequences in a fasta file.

4. `count_bases_whole_genome.py` -- the script counts all IUPAC bases in genome sequences in .gbk files

5. `merge_bases_categories_taxonomy` -- the script merges all info about gene sequences (nucleotide composition, category, taxonomy) into a single TSV file.

6. `seq_logos/logo_making.ipynb` -- a Jupyter notebook used for sequence logo creation.

7. `alt_antiSD_tests.R` -- a script used for testing is a "mixed" 16S rRNA gene set is uniformly distributed accross the bacterial phyla.


