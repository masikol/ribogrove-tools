# RiboGrove paper

This repo contains the software for RiboGrove database constuction and exploration.

RiboGrove is a genome-based database of full-length sequences of prokaryotic 16S rRNA genes. Sequences are derived from RefSeq database.

## Latest release

(not yet) The first RiboGrove release (1.207) was constructed using scripts from this repo released under version 1.207 (as well). RefSeq release 207 was used for that, hence the RiboGrove release version.

## Contents:

1. Scripts in direcory `db_creation_and_filtering/` are used for collection of 16S rRNA genes from RefSeq, and for further removal (i.e. filtering) of partial gene sequencing.

2. Scripts in direcory `exploration_scripts/` are used for primer coverage calculation, calculation of intragenomic variability of 16S rRNA genes etc.

3. Demonstration data in the directory `demo/`. This is how input data and a working directory for the filtering pipeline looks like.

4. (scripts in directory `_trash/` are deprecated)

## Python version

All Python scripts in this repo are written for Python 3 (version 3.6 or later).

