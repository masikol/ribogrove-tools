
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

## make_dada2_train_set_file.py

### Input

This script takes RiboGrove fasta file `ribogrove_XX.XXX_sequences.fasta.gz` as input (`XX.XXX` is RiboGrove release number, e.g. 10.216).

### Usage:

Without specific epithets:

```
python3 make_dada2_train_set_file.py \
    -i ribogrove_10.216_sequences.fasta.gz \
    -o ribogrove_10.216_dada2_train_set.fasta.gz
```

With specific epithets (set `-w` option):

```
python3 make_dada2_train_set_file.py \
    -i ribogrove_10.216_sequences.fasta.gz \
    -o ribogrove_10.216_dada2_train_set.fasta.gz \
    -w
```
