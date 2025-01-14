
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
  --base-counts      <BASE_COUNTS_TSV> \
  --taxonomy         <TAXONOMY_TSV> \
  --categories       <CATEGORIES_TSV> \
  --entropy-summary  <ENTROPY_SUMMARY_TSV> \
  --source-genomes   <SOURCE_GENOMES_TSV>  \
  --primers-dir      <PRIMER_COVERAGE_DIR> \
  --outdir           <OUTDIR>              \
  --seqkit           <SEQKIT_EXECUTABLE>
```

### Example

```
python3 make_ribogrove_release_page.py \
  --release-num 22.228 \
  --release-date 2025-01-14 \
  --final-fasta /mnt/cager-beast/m.2.2/RiboGrove/RiboGrove_releases/22.228/ribogrove_22.228_sequences.fasta.gz \
  --metadata /mnt/cager-beast/m.2.2/RiboGrove/RiboGrove_releases/22.228/metadata_ribogrove_22.228.zip \
  --base-counts /mnt/cager-beast/m.2.2/RiboGrove/RiboGrove_releases/22.228/metadata/gene_seqs_base_counts.tsv \
  --taxonomy /mnt/cager-beast/m.2.2/RiboGrove/RiboGrove_releases/22.228/metadata/taxonomy.tsv \
  --categories /mnt/cager-beast/m.2.2/RiboGrove/RiboGrove_releases/22.228/metadata/categories.tsv \
  --entropy-summary /mnt/cager-beast/m.2.2/RiboGrove/RiboGrove_releases/22.228/metadata/entropy_summary.tsv \
  --source-genomes /mnt/cager-beast/m.2.2/RiboGrove/RiboGrove_releases/22.228/metadata/source_RefSeq_genomes.tsv \
  --primers-cov /mnt/cager-beast/m.2.2/RiboGrove/RiboGrove_releases/22.228/metadata/primer_pair_genomic_coverage.tsv \
  --outdir /mnt/cager-beast/m.2.2/RiboGrove/RiboGrove_pages/22.228 \
  --seqkit /usr/bin/seqkit
```


## Dependencies

### Python packages

1. numpy. Installation: `pip3 install numpy`. Tested on version 2.2.1.

2. pandas. Installation: `pip3 install pandas`. Tested on version 2.2.3.

3. flask. Installation: `pip3 install Flask`. Tested on version 3.1.0.

### Standalone programs

1. [seqkit](https://github.com/shenwei356/seqkit). Tested on version 2.9.0.
