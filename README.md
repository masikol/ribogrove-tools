# RiboGrove tools

This repo contains the software for RiboGrove database constuction and exploration.

[RiboGrove](http://mbio.bas-net.by/cager/en/ribogrove) is a genome-based database of full-length sequences of prokaryotic 16S rRNA genes. Sequences are derived from RefSeq database.

The [main website](http://mbio.bas-net.by/cager/en/ribogrove) where RiboGrove is hosted may be unavailable outside Belarus due to technical troubles and the overall disaster. Nevertheless, RiboGrove releases are available at a Github Pages mirror: [https://masikol.github.io/](https://masikol.github.io/). We are sorry for this inconvenience.

## Contents:

1. Scripts in directory `db_creation_and_filtering/` are used for collection of 16S rRNA genes from RefSeq, and for further removal (i.e. filtering) of partial gene sequencing.

2. Scripts in directory `release_preparation/` are used for making release files (`*.fasta.gz` and `metadata.zip`) and HTML release pages.

3. Scripts in directory `exploration_scripts/` are used for draw some plots and calculate some additional data for the publication etc.

4. Demonstration data in the directory `demo/`. This is how input data and a working directory for the filtering pipeline looks like.

5. (scripts in directory `_trash/` are deprecated)

## Releases of ribogrove-tools

Here is a list of RiboGrove database releases and corresponding ribogrove-tools releases the database was created with.

- RiboGrove database releases 1.207–7.213 -- ribogrove-tools release `1.207`.

- RiboGrove database release 8.214 -- ribogrove-tools release `8.214`.

- RiboGrove database release 9.215 -- ribogrove-tools release `9.215`.

## Python version

All Python scripts in this repo are written for Python 3 (version 3.6 or later).

## Citing RiboGrove

If you find RiboGrove useful for your research please cite:

Maxim A. Sikolenko, Leonid N. Valentovich. “RiboGrove: a database of full-length prokaryotic 16S rRNA genes derived from completely assembled genomes” // Research in Microbiology, Volume 173, Issue 4, May 2022, 103936.
(DOI: [10.1016/j.resmic.2022.103936](https://doi.org/10.1016/j.resmic.2022.103936))
