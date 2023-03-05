# Scripts for RiboGrove creation pipeline

Here you can find scripts implementing individual steps of RiboGrove creation pipeline:

## Download annotated genome sequences in GenBank format from RefSeq

1. `filter_refseq_catalog.py` filters `.catalog` file of the target RefSeq release. “Filters” means removes records which don’t belong to [“Genomic” molecule type](https://www.ncbi.nlm.nih.gov/books/NBK21091/table/ch18.T.refseq_accession_numbers_and_mole/?report=objectonly/).

2. Then, [Assembly Summary](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/${DOMAIN}/assembly_summary.txt) file gets downloaded.

3. `filter_asm_summary_step1.py` filters the downloaded Assembly Summary: 1) retains only genomes of “Complete Genome” and “Chromosome” completeness; 2) removes “Whole genome shotgun” sequences, since they don’t represent completely assembled genomes; 3) removes genomes from the blacklist (`data/ad_hoc/assembly_blacklist.tsv`).

4. `download_genomes.py` downloads genomes: 1) `.gbff.gz` files of annotated genomic sequences; 2) `_assembly_report.txt` files with brief info such as sequencing technologies.

5. `make_replicon_map.py` makes a helper file in which Assembly acceession numbers (for example `GCF_000005825.2`) are mapped to corresponding RefSeq accession numbers (in the example, `NC_013791.2`, `NC_013792.1` and `NC_013793.1`).

6. `filter_asm_summary_step2.py`, having genome sequences and RefSeq accessions, filters downloaded Assembly Summary again: 1) removes genomes which don’t belong to the current RefSeq release using the `.catalog` file; 2) removes genomes with sequences containing at least 3 Ns in row.

7. `make_taxonomy.py` takes a `new_taxdump/rankedlineage.dmp` file and assigns taxonomy to the downloaded and filtered genomes.

## Extract 16S rRNA gene sequences from downloaded genomes

8. `extract_16S.py` extracts sequences of 16S rRNA genes from the downloaded and filtered genomes.

## Assign categories to the genomes according to their reliability

9. `assign_genome_categories.py` assigns categories to the genomes according to their reliability: from 1 (for high reliability) to 3 (for low reliability).

## Filter gene sequences: remove partial sequences, sequences containing long repeats etc.

10. `check_seqs_with_ribotyper.py` processes 16S rRNA gene sequences with `ribotyper` to find unreliable sequences.

11. `find_ribotyper_fail_seqs.py` records which gene sequences didn’t pass the `ribotyper` filter, i.e. sequences with the following unexpected features: “NoHits”, “UnacceptableModel”, “MinusStrand”, “LowScore”, “LowCoverage”.

12. `find_aberrant_genes.py` finds partial genes and genes having large deletion compared to the corresponding pivotal gene. A pivotal gene is a gene having maximum `tscore` determined by `ribotyper` among all 16S rRNA genes from a genome.

13. `find_repeats.py` finds gene sequences having long internal repeats.

14. `make_final_seqs.py` removes gene sequences which failed at least one of the three filters: “find_ribotyper_fail_seqs”, “find_aberrant_genes” and “find_repeats”.

## Annotate final gene sequences and prepare additional files

15. `annotate_seq_names.py` adds information to headers of final sequences, namely 1) taxonomy; 2) genome category.

16. `count_bases.py` tallies bases (i.e. letters `AGCT` and degenerate `RYWSMKHVBDN`) in final sequences.

17. `merge_bases_categories_taxonomy` aggregates 1) taxonomy; 2) genome categories; 3) bases counts for each gene sequence into a single file.

18. `calculate_entropy.py` calculates intragenomic variability of final gene sequences.

19. `check_primers_mfeprimer.py` calculates coverage of different primer pairs to different V-regions of bacterial 16S rRNA genes.

20. `calculate_GCNs.py` calculates 16S rRNA Gene Copy Numbers per genome. Also, it can calculate "primer-wise" GCNs. In this mode, the script also tests if PCR primers anneal to gene sequences. Thus, the script will give a genome +1 copy only if the gene sequence can form a PCR product with a primer pair.