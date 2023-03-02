
# Data transformation scripts

These scripts can help you convert RiboGrove release files to several external formats:

1. To QIIME2-compatible taxonomy file (script `make_qiime_taxonomy_file.py`).

2. To a reference fasta file for dada2 taxonomic classificator (script `make_dada2_train_set_file.py`).


## make_qiime_taxonomy_file.py

### Input

This script takes file `gene_seqs_statistics.tsv` from RiboGrove metadata as input.

### Usage:

```
python3 make_qiime_taxonomy_file.py \
    -i ribogrove_release/metadata/gene_seqs_statistics.tsv \
    -o qiime2_compatible_taxonomy.tsv
```

### Dependencies

#### Python packages

1. numpy. Installation: `pip3 install numpy`. Tested on version TODO: add version.

2. pandas. Installation: `pip3 install pandas`. Tested on version 1.1.3.

## make_dada2_train_set_file.py

### Input

This script takes RiboGrove fasta file `ribogrove_XX.XXX_sequences.fasta.gz` as input (`XX.XXX` is RiboGrove release number, e.g. 10.216).

### Usage:

Make a dada2-train-set file and include specific epithets in taxonomy strings:

```
python3 make_dada2_train_set_file.py \
    -i ribogrove_10.216_sequences.fasta.gz \
    -o ribogrove_10.216_dada2_train_set.fasta.gz \
    -w
```

Without specific epithets -- simply don’t set `-w` option:

```
python3 make_dada2_train_set_file.py \
    -i ribogrove_10.216_sequences.fasta.gz \
    -o ribogrove_10.216_dada2_train_set.fasta.gz \
```

### Dependencies

#### Python packages

1. Biopython. Installation: `pip3 install biopython`. Tested on version TODO: add version.
