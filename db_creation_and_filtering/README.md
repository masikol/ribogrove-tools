# RiboGrove database creation: sequence collection and filtering

Here you can find scripts that are used for RiboGrove construction.

Each Python script can be used on its own using command line interace. At the same time, they are integrated in a pipeline `collect_and_filter_ssu_pipeline.sh`.

## Description

The pipeline consists of 5 major steps:

1. Download annotated genome sequences in GenBank format from RefSeq.

2. Extract 16S rRNA gene sequences from downloaded genomes.

3. Collect NCBI taxonomy data for the downloaded genomes.

4. Filter gene sequences: remove partial sequences.

5. Annotate final gene sequences and prepare additional files.

The pipeline must be run sparately for bacteria and for archaea.

## Input data for the pipeline

1. A text file listing Assembly IDs from NCBI [Assembly](https://www.ncbi.nlm.nih.gov/assembly/) database, one per line. Genomes sequences of listed assemblies are used for RiboGrove creation.

2. A RefSeq `.catalog.gz` file from [https://ftp.ncbi.nlm.nih.gov/refseq/release/release-catalog/](https://ftp.ncbi.nlm.nih.gov/refseq/release/release-catalog/). E.g. for RiboGrove 7.213 you should use file `RefSeq-release213.catalog.gz`.

3. A `rankedlineage.dmp` file from NCBI [Taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy) database. This file can be downloaded within `new_taxdump` archive from [NCBI Taxonomy FTP site](https://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/). It is better to download this archive before construction of each RiboGrove release to use the most recent taxonomic data.

The pipeline takes a bash-like configuration file (`.conf`) as input. Examples of that configuration files can be found in directory `config/`: `bacteria_example.conf` and `archaea_example.conf`.

### Usage

Usage of the pipeline script:

```bash
bash collect_and_filter_ssu_pipeline.sh <CONFIG_FILE>
```

Example:

```bash
bash collect_and_filter_ssu_pipeline.sh config/bacteria_config.conf
```

### Configuragion file

A configuragion file should specify locations of 1) input data; 2) working directory; 3) dependencies. Also, the configuration file specifies some pipeline parameters, such as Rfam RNA family ID. Bash environment variables are used for this purpose.

#### 1. Input data, output location

1. `PREFIX` -- prefix for output files (e.g. `archaea`, `bacteria`).

2. `ASSEMBLY_IDS_FPATH` -- a file listing [Assembly](https://www.ncbi.nlm.nih.gov/assembly/) IDs, one per line.

3. `GENOMES_GBK_DIR` -- directory for `.gbk.gz` files of downloaded genomes.

4. `RANKEDLINEAGE_FPATH` -- file `rankedlineage.dmp` from archive `new_taxdump`. This file contains NCBI taxonomy.

#### 2. Working directory

1. `WORKDIR` -- working directory, i.e. directory for output and temporary files.

#### 3. Dependencies

1. `SEQKIT` -- [seqkit](https://github.com/shenwei356/seqkit) executable.

2. `REFSEQ_CATALOG_FILE` -- a A RefSeq `.catalog.gz` file described [above](#input-data-for-the-pipeline).

3. `CMSEARCH_FOR_EXTRACT_16S` -- the program `cmsearch` from [Infernal](http://eddylab.org/infernal/). The program is used to extract 16S rRNA genes from genomic sequences. Here, we use `cmsearch` 1.1.1, to be fully consistent wuth PGAP.

4. `RFAM_FOR_EXTRACT_16S` -- [Rfam](https://rfam.xfam.org/) covariance model (uncompressed `.cm` file). It is used to extract 16S rRNA genes from genomic sequences. Here, we use Rfam 12.0: ftp://ftp.ebi.ac.uk/pub/databases/Rfam/12.0/Rfam.cm.gz, to be fully consistent wuth PGAP.

5. `CMSCAN_FOR_FILTERING` -- the program `cmscan` from [Infernal](http://eddylab.org/infernal/). The program is used to filter incomplete sequences of 16S rRNA genes. Here, we use the latest version of `cmscan` (1.1.4).

6. `CMPRESS_FOR_FILTERING` -- the program `cmpress` from [Infernal](http://eddylab.org/infernal/). The program is used to filter incomplete sequences of 16S rRNA genes. Here, we use the latest version of `cmpress` (1.1.4).

7. `RFAM_FOR_FILTERING` -- [Rfam](https://rfam.xfam.org/) covariance model (uncompressed `.cm` file). It is used to extract 16S rRNA genes from genomic sequences. Here, we use latest version Rfam 14.6: ftp://ftp.ebi.ac.uk/pub/databases/Rfam/14.6/Rfam.cm.gz.

8. `CMFETCH` -- the program `cmfetch` from [Infernal](http://eddylab.org/infernal/). The program is used to extract specific Rfam model from from whole Rfam file. Here, we use the latest version of `cmfetch` (1.1.4).

9. `RFAM_FAMILY_ID` -- ID of specific Rfam covariance model to use. We used the following IDs to predict 16S rRNA: `RF00177` for bacteria and `RF01959` for archaea.

10. `MUSCLE` -- [MUSCLE](https://www.drive5.com/muscle/) executable for multiple sequence alignment.

1. `MFEPRIMER` -- [MFEprimer](https://www.mfeprimer.com/) executable for primer coverage calculations.

#### 4. Conserved 16S rRNA regions

1. `CONSERVED_REGIONS_FASTA` -- a fasta file of the 16S rDNA conserved regions. We use the bundled file `db_creation_and_filtering/conserved_regions.fasta`. These sequences are adopted from [M. Martinez-Porchas et al](https://doi.org/10.7717/peerj.3036).

2. `CHECK_CONSERV_REGIONS` -- (accepted values: `0`, `1`). `1` -- check the 16S rDNA conserved regions while finding aberrant genes. `0` -- do not check the conserved regions; detect aberrant genes only by deletions instead.

#### 5. Primer coverage calculation

1. `CALC_PRIMERS_COVERAGE` -- (accepted values: `0`, `1`). `1` -- calculate primer coverage with MFEprimer. `0` -- do not calculate primer coverage.

#### 6. "Cached" data from previous RiboGrove release

1. `PREV_WORKDIR` -- a working directory of the previous RiboGrove release. If you are running constructing RiboGrove 8.214 for bacteria, your `PREV_WORKDIR` value will be something like `/some/path/RiboGrove_workdirs/7.213/bacteria`.

## Dependencies

Python packages:

1. numpy. Installation: `pip3 install numpy`. Tested on version 1.19.2.

2. pandas. Installation: `pip3 install pandas`. Tested on version 1.1.3.

3. Biopython. Installation: `pip3 install biopython`. Tested on version 1.78.

4. repeatfinder ([repo](https://github.com/deprekate/RepeatFinder)). Installation: `pip3 install repeatfinder`. Tested on version 1.5.

Standalone programs:

1. Seqkit: [https://github.com/shenwei356/seqkit](https://github.com/shenwei356/seqkit). Tested on version 2.2.0.

2. MUSCLE: [https://www.drive5.com/muscle/](https://www.drive5.com/muscle/). Tested on version 3.8.31.

3. Infernal [http://eddylab.org/infernal/](http://eddylab.org/infernal/). Infernal 1.1.1 is used for re-annotation of genomic sequences (for compatibility with [PGAP](https://www.ncbi.nlm.nih.gov/genome/annotation_prok/)), but in further steps newer versions of Infernal may be used for sequence filtering (e.g. version 1.1.4).

3. MFEprimer [https://www.mfeprimer.com/](https://www.mfeprimer.com/). We use MFEprimer version 3.2.2.

Data:

1. Rfam database. Versions: 12.0 (ftp://ftp.ebi.ac.uk/pub/databases/Rfam/12.0/Rfam.cm.gz) for re-annotation of genomic sequences; 14.6 (ftp://ftp.ebi.ac.uk/pub/databases/Rfam/14.6/Rfam.cm.gz) for sequence filtering.
