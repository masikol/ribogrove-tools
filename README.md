# RiboGrove tools

This repo contains the software for RiboGrove database constuction and exploration.

[RiboGrove](mbio.bas-net.by/cager/en/ribogrove) is a genome-based database of full-length sequences of prokaryotic 16S rRNA genes. Sequences are derived from RefSeq database.

The [main website](mbio.bas-net.by/cager/en/ribogrove) where RiboGrove is hosted may be unavailable outside Belarus due to technical troubles. We will host the latest RiboGrove release at Github Pages: [https://masikol.github.io/](https://masikol.github.io/). We are sorry for this inconvenience.

## Releases of ribogrove-tools

The first [RiboGrove release](http://mbio.bas-net.by/cager/en/content/59-ribogrove-1-207) (1.207) was constructed using scripts from this repo released under version 1.207. RefSeq release 207 was used for that, hence the RiboGrove release version.

The software at the `main` branch are **not** for public use, only the releases are.

## Contents:

1. Scripts in directory `db_creation_and_filtering/` are used for collection of 16S rRNA genes from RefSeq, and for further removal (i.e. filtering) of partial gene sequencing.

2. Scripts in directory `release_preparation/` are used for making release files (`.fasta.gz` and `metadata.zip`), and HTML release pages.

3. Scripts in directory `exploration_scripts/` are used for primer coverage calculation, calculation of intragenomic variability of 16S rRNA genes etc.

4. Demonstration data in the directory `demo/`. This is how input data and a working directory for the filtering pipeline looks like.

5. (scripts in directory `_trash/` are deprecated)

## Python version

All Python scripts in this repo are written for Python 3 (version 3.6 or later).

## Citing RiboGrove

If you find RiboGrove useful for your research please cite:

Maxim A. Sikolenko, Leonid N. Valentovich. "RiboGrove: a database of full-length prokaryotic 16S rRNA genes derived from completely assembled genomes" // Research in Microbiology, 2022, 103936.
(DOI: [10.1016/j.resmic.2022.103936](https://doi.org/10.1016/j.resmic.2022.103936). Epub ahead of print).
