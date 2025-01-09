# RiboGrove database creation: sequence collection and filtering

Here you can find scripts that are used for RiboGrove construction.

Each Python script can be used on its own using command line interface. At the same time, they are integrated into a pipeline `collect_and_filter.sh`.

## Description

The pipeline consists of 5 major steps:

1. Download annotated genome sequences in GenBank format from RefSeq.

2. Extract 16S rRNA gene sequences from downloaded genomes.

3. Assign categories to the genomes according to their reliability.

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

The pipeline should be run as follows:

```
bash collect_and_filter.sh <CONFIG_FILE> [-at]
```

### Options

```
-t: Test mode. If specified, the script will not use entire RefSeq,
    but just some genomes from the file
    scripts/data/test/test_bacteria_assembly_summary.txt.gz
    or
    scripts/data/test/test_archaea_assembly_summary.txt.gz

-a: RefSeq catalog file is already filtered.
    If you already have file RefSeq-release216_filtered.catalog.gz
    created from file RefSeq-release216.catalog.gz,
    the script will not run the script filter_refseq_catalog.py again
    and thus save some time.
```

### Example:

```
bash collect_and_filter.sh config/bacteria_config.conf -t
```

## Dependencies

### Python packages:

1. numpy. Installation: `pip3 install numpy`. Tested on version 1.19.2.

2. pandas. Installation: `pip3 install pandas`. Tested on version 1.1.3.

3. Biopython. Installation: `pip3 install biopython`. Tested on version 1.81.

4. [repeatfinder](https://github.com/deprekate/RepeatFinder). Installation: `pip3 install repeatfinder`. Tested on version 1.5.

### Standalone programs:

1. Seqkit: [https://github.com/shenwei356/seqkit](https://github.com/shenwei356/seqkit). Tested on version 2.9.0.

2. `ribotyper` program from [Ribovore](https://github.com/ncbi/ribovore) suite. Tested on Ribovore version 1.0.5.

2. MAFFT: [https://mafft.cbrc.jp/alignment/software](https://mafft.cbrc.jp/alignment/software). Tested on version v7.526.

3. Infernal [http://eddylab.org/infernal/](http://eddylab.org/infernal/). Infernal 1.1.5 is used for re-annotation of genomic sequences (for compatibility with [PGAP](https://www.ncbi.nlm.nih.gov/genome/annotation_prok/)).

3. MFEprimer [https://www.mfeprimer.com/](https://www.mfeprimer.com/). Tested on version 3.3.1.

### Data:

1. Rfam database version 14.4 (ftp://ftp.ebi.ac.uk/pub/databases/Rfam/14.4/Rfam.cm.gz) for re-annotation of genomic sequences for compatibility with [PGAP](https://www.ncbi.nlm.nih.gov/genome/annotation_prok/).
