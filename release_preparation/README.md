# RiboGrove release preparation

The software in this directory should be used to create a RiboGrove release. If you have output directories made by the pipeline `collect_and_filter_ssu_pipeline.sh`, you can the create a release using scripts in this directory.

If the release is created, you can then make HTML pages in different languages for this release.

## Release preparation

You should use script `prepare_release.sh` for RiboGrove release preparation. It will crate two fasta files: a full-length RiboGrove file and a "raw" version of it, i.e. without filtering. Also, the script will create a metadata archive with varuios additional information.

### Usage

```bash
bash prepare_release.sh \
    <WORKDIR> \
    <OUTDIR> \
    <RELEASE_NUMBER>
```

#### Example

```bash
bash prepare_release.sh \
    /mnt/1.5_drive_0/RiboGrove/RiboGrove_workdirs/7.213 \
    /mnt/1.5_drive_0/RiboGrove/RiboGrove_releases/7.213 \
    7.213
```

`<WORKDIR>` directory should be of the following structure:

```
├── archaea/
└── bacteria/
```

`archaea/` and `bacteria/` directories are working directories of the pipeline `collect_and_filter_ssu_pipeline.sh`. You set their names using `PREFIX` environment variable of configuration file for the pipeline. Thus, `prepare_release.sh` expects these prefixes to be `archaea` and `bacteria`.

### Dependencies

This script depends on two programs:

1. [csvtk](https://github.com/shenwei356/csvtk). Tested on version 0.24.0.

2. [seqkit](https://github.com/shenwei356/seqkit). Ttested on version 2.2.0.

The programs should be in PATH variable.


## Release page preparation

You should use script `make_ribogrove_release_page.py` to create HTML release pages for RiboGrove.

### Usage:

```bash
python3 make_ribogrove_release_page.py \
  --release-num      <RELEASE_NUMBER>      \
  --release-date     <RELEASE_DATE>        \
  --final-fasta      <FULL_LENGTH_FASTA>   \
  --raw-fasta        <RAW_LENGTH_FASTA>    \
  --metadata         <METADATA_ZIP>        \
  --gene-stats-table <GENE_STATS_TSV>      \
  --entropy-summary  <ENTROPY_SUMMARY_TSV> \
  --source-genomes   <SOURCE_GENOMES_TSV>  \
  --outdir           <OUTDIR>              \
  --seqkit           <SEQKIT_EXECUTABLE>
```

#### Example

```bash
python3 make_ribogrove_release_page.py \
  --release-num      7.213 \
  --release-date     2022-07-20 \
  --final-fasta      /mnt/1.5_drive_0/RiboGrove_releases/7.213/ribogrove_7.213_sequences.fasta.gz \
  --raw-fasta        /mnt/1.5_drive_0/RiboGrove_releases/7.213/raw_ribogrove_7.213_sequences.fasta.gz \
  --metadata         /mnt/1.5_drive_0/RiboGrove_releases/7.213/metadata_ribogrove_7.213.zip \
  --gene-stats-table /mnt/1.5_drive_0/RiboGrove_releases/7.213/metadata/gene_seqs_statistics.tsv \
  --entropy-summary  /mnt/1.5_drive_0/RiboGrove_releases/7.213/metadata/entropy_summary.tsv \
  --source-genomes   /mnt/1.5_drive_0/RiboGrove_releases/7.213/metadata/source_RefSeq_genomes.tsv \
  --outdir           /mnt/1.5_drive_0/RiboGrove_pages/7.213/ \
  --seqkit           /usr/bin/seqkit
```

### Dependencies

#### Python packages

1. pandas. Installation: `pip3 install pandas`. Tested on version 1.1.3.

2. flask. Installation: `pip3 install Flask`. Tested on version 2.0.2.

#### Standalone programs

1. [seqkit](https://github.com/shenwei356/seqkit). Tested on version 2.2.0.