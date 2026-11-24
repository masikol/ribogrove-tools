
# RiboGrove HTML pages

The script `make_ribogrove_release_page.py` creates HTML pages for RiboGrove releases.

## Input

RiboGrove release, namely files `ribogrove_XX.XXX_sequences.fasta.gz` and `metadata.zip` (`XX.XXX` is RiboGrove release number, e.g. 10.216).

Also, the script requires directory for primer coverages of bacterial genomes (option `--primers-dir`).

## Usage:

```
python3 make_ribogrove_release_page.py \
  --release-num      <RELEASE_NUMBER>       \
  --release-date     <RELEASE_DATE>         \
  --final-fasta      <FULL_LENGTH_FASTA>    \
  --metadata         <METADATA_ZIP>         \
  --base-counts      <BASE_COUNTS_TSV>      \
  --taxonomy         <TAXONOMY_TSV>         \
  --categories       <CATEGORIES_TSV>       \
  --entropy-summary  <ENTROPY_SUMMARY_TSV>  \
  --source-genomes   <SOURCE_GENOMES_TSV>   \
  --primers-cov      <PRIMERS_COV_TSV>      \
  --density          <16S_GENE_DENSITY_TSV> \
  --zenodo-doi       <ZENODO_DOI>           \
  --zenodo-record-id <ZENODO_RECORD_ID>     \
  --outdir           <OUTDIR>               \
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
  --density /mnt/cager-beast/m.2.2/RiboGrove/RiboGrove_releases/22.228/metadata/16S_density.tsv \
  --zenodo-doi 10.5281/zenodo.17200498 \
  --zenodo-record-id 17200498 \
  --outdir /mnt/cager-beast/m.2.2/RiboGrove/RiboGrove_pages/22.228 \
  --seqkit /usr/bin/seqkit
```


## Dependencies

### Python packages

1. numpy. Installation: `pip3 install numpy`. Tested on version 2.5.3.

2. polars. Installation: `pip3 install polars`. Tested on version 2.2.3.

3. flask. Installation: `pip3 install Flask`. Tested on version 3.1.0.

### Standalone programs

1. [seqkit](https://github.com/shenwei356/seqkit). Tested on version 2.9.0.
