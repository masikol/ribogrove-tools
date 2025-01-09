# Prepare RiboGrove release

The script `combine_release.sh` (found this directory) combines bacterial and archaeal data into a RiboGrove release, namely into files `ribogrove_XX.XXX_sequences.fasta.gz` and `metadata.zip` (`XX.XXX` is RiboGrove release number, e.g. 10.216).

## Input

Output working directories `bacteria/` and `archaea/` made by the pipeline `collect_and_filter/collect_and_filter.sh`.

## Usage

```
bash combine_release.sh \
    <WORKDIR_ROOT> \
    <OUTDIR> \
    <RELEASE_NUMBER>
```

### Example

```
bash combine_release.sh \
    /mnt/my_drive/RiboGrove_workdirs/9.215 \
    /mnt/my_drive/RiboGrove_releases/9.215 \
    9.215
```

`<WORKDIR_ROOT>` directory should be of the following structure:

```
├── archaea/
└── bacteria/
```

`archaea/` and `bacteria/` directories are working directories of the pipeline `collect_and_filter/collect_and_filter.sh`.

### Dependencies

This script depends on two programs:

1. [csvtk](https://github.com/shenwei356/csvtk). Tested on version 0.24.0.

2. [seqkit](https://github.com/shenwei356/seqkit). Tested on version 2.2.0.

The programs should be in PATH environment variable.
