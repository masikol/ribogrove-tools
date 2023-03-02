# RiboGrove database creation: sequence collection and filtering

Here you can find scripts that are used for RiboGrove construction.

Each Python script can be used on its own using command line interface. At the same time, they are integrated into a pipeline `collect_and_filter.sh`.

## Description

The pipeline consists of 5 major steps:

1. Download annotated genome sequences in GenBank format from RefSeq.

2. Extract 16S rRNA gene sequences from downloaded genomes.

3. Assign categories to the genomes accorging to their reliability.

4. Filter gene sequences: remove partial sequences, sequences containing long repeats etc.

5. Annotate final gene sequences and prepare additional files.

The pipeline must be run sparately for bacteria and for archaea.

## Input

## Configuragion file

The pipeline takes a bash-like configuration file (`.conf`) as input. Examples of that configuration files can be found in directory `config/`: `bacteria_example.conf` and `archaea_example.conf`.

A configuration file specifies location of:

1) input data;

2) working directory;

3) dependencies.

Please see examples of configuration files in the directory `config/`.

Also, the configuration file specifies some pipeline parameters, such as Rfam RNA family ID and whether to calculate primer coverages of genomes or not.

## NCBI files to download before run

The pipeline requires the following files downloaded:

1. A RefSeq `.catalog.gz` file from [https://ftp.ncbi.nlm.nih.gov/refseq/release/release-catalog/](https://ftp.ncbi.nlm.nih.gov/refseq/release/release-catalog/). E.g. for RiboGrove 7.213 you should use file `RefSeq-release213.catalog.gz`.

2. A `rankedlineage.dmp` file from NCBI [Taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy) database. This file can be downloaded within `new_taxdump` archive from [NCBI Taxonomy FTP site](https://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/). It is better to download this archive before construction of each RiboGrove release to use up-to-date taxonomic data.

## Usage

Usage of the pipeline script:

```
bash collect_and_filter.sh <CONFIG_FILE> [-at]
```

Example:

```
bash collect_and_filter.sh config/bacteria_config.conf -t
```

### Options

```
-t: TODO: add help

-a: TODO: add help
```

## Dependencies

### Python packages:

1. numpy. Installation: `pip3 install numpy`. Tested on version 1.19.2.

2. pandas. Installation: `pip3 install pandas`. Tested on version 1.1.3.

3. Biopython. Installation: `pip3 install biopython`. Tested on version 1.78.

4. repeatfinder ([repo](https://github.com/deprekate/RepeatFinder)). Installation: `pip3 install repeatfinder`. Tested on version 1.5.

### Standalone programs:

1. Seqkit: [https://github.com/shenwei356/seqkit](https://github.com/shenwei356/seqkit). Tested on version 2.2.0.

2. `ribotyper` program from [Ribovore](https://github.com/ncbi/ribovore) suite. Tested on version TODO: add version.

2. MUSCLE: [https://www.drive5.com/muscle/](https://www.drive5.com/muscle/). Tested on version 3.8.31.

3. Infernal [http://eddylab.org/infernal/](http://eddylab.org/infernal/). Infernal 1.1.1 is used for re-annotation of genomic sequences (for compatibility with [PGAP](https://www.ncbi.nlm.nih.gov/genome/annotation_prok/)).

3. MFEprimer [https://www.mfeprimer.com/](https://www.mfeprimer.com/). Tested on version 3.2.2.

### Data:

1. Rfam database. Versions: 12.0 (ftp://ftp.ebi.ac.uk/pub/databases/Rfam/12.0/Rfam.cm.gz) for re-annotation of genomic sequences.
