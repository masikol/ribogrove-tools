
# RiboGrove HTML pages

The script `make_ribogrove_release_page.py` creates HTML pages for RiboGrove releases.

## Input

RiboGrove release, namely files `ribogrove_XX.XXX_sequences.fasta.gz` and `metadata.zip` (`XX.XXX` is RiboGrove release number, e.g. 10.216).

Also, the script requires directory for primer coverages of bacterial genomes (option `--primers-dir`).

## Usage:

```
python3 make_ribogrove_release_page.py \
  --release-num      <RELEASE_NUMBER>      \
  --release-date     <RELEASE_DATE>        \
  --final-fasta      <FULL_LENGTH_FASTA>   \
  --metadata         <METADATA_ZIP>        \
  --gene-stats-table <GENE_STATS_TSV>      \
  --entropy-summary  <ENTROPY_SUMMARY_TSV> \
  --source-genomes   <SOURCE_GENOMES_TSV>  \
  --primers-dir      <PRIMER_COVERAGE_DIR> \
  --outdir           <OUTDIR>              \
  --seqkit           <SEQKIT_EXECUTABLE>
```

### Example

```
python3 make_ribogrove_release_page.py \
  --release-num      9.215 \
  --release-date     2022-07-20 \
  --final-fasta      /mnt/my_drive/RiboGrove_releases/9.215/ribogrove_9.215_sequences.fasta.gz \
  --metadata         /mnt/my_drive/RiboGrove_releases/9.215/metadata_ribogrove_9.215.zip \
  --gene-stats-table /mnt/my_drive/RiboGrove_releases/9.215/metadata/gene_seqs_statistics.tsv \
  --entropy-summary  /mnt/my_drive/RiboGrove_releases/9.215/metadata/entropy_summary.tsv \
  --source-genomes   /mnt/my_drive/RiboGrove_releases/9.215/metadata/source_RefSeq_genomes.tsv \
  --primers-dir      /mnt/my_drive/RiboGrove/RiboGrove_workdirs/9.215/bacteria/primers_coverage \
  --outdir           /mnt/my_drive/RiboGrove_pages/9.215/ \
  --seqkit           /usr/bin/seqkit
```


## Dependencies

### Python packages

1. numpy. Installation: `pip3 install numpy`. Tested on version 1.19.2.

2. pandas. Installation: `pip3 install pandas`. Tested on version 1.1.3.

3. flask. Installation: `pip3 install Flask`. Tested on version 2.0.2.

### Standalone programs

1. [seqkit](https://github.com/shenwei356/seqkit). Tested on version 2.9.0.
