# RiboGrove paper

This repo contains the software for RiboGrove database constuction and exploration.

RiboGrove is a genome-based database of full-length sequences of prokaryotic 16S rRNA genes. Sequences are derived from RefSeq database.

## Releases of ribogrove-paper

The first [RiboGrove release](http://mbio.bas-net.by/cager/en/content/59-ribogrove-1-207) (1.207) was constructed using scripts from this repo released under version 1.207. RefSeq release 207 was used for that, hence the RiboGrove release version.

The ribogrove-paper releases for RiboGrove 2.208 and 3.209 are coming soon.

The software at the `main` branch are **not** for public use, only the releases are.

## Contents:

1. Scripts in directory `db_creation_and_filtering/` are used for collection of 16S rRNA genes from RefSeq, and for further removal (i.e. filtering) of partial gene sequencing.

2. Scripts in directory `release_preparation/` are used for making release files (`.fasta.gz` and `metadata.zip`), and HTML release pages.

3. Scripts in directory `exploration_scripts/` are used for primer coverage calculation, calculation of intragenomic variability of 16S rRNA genes etc.

4. Demonstration data in the directory `demo/`. This is how input data and a working directory for the filtering pipeline looks like.

5. (scripts in directory `_trash/` are deprecated)

## Python version

All Python scripts in this repo are written for Python 3 (version 3.6 or later).

