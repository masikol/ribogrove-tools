# RybaSom database creation and sequence filtering

Here you can find scripts that were (and can be) used for RybaSom construction.

Each Python script can be used on its own using command line interace. At the same time, they are integrated in a pipeline `collect_and_filter_ssu_pipeline.sh`.

## Pipeline description and usage

### Description

The pipeline consists of 4 major steps:

1. Download annotated genome sequences in GenBank format from RefSeq.

2. Extract 16S rRNA gene sequences from downloaded genomes.

3. Remove partial sequences (i.e. perform filtering).

4. Annotate result sequences with taxonomy data and genome categories.

### Input data for the pipeline

The pipeline takes configuration file (`.conf`) as input. Examples of that configuration files can be found in directory `config/`: `bacteria_config.conf` and `archaea_config.conf`.

Thus, usage of pipeline script is following:

```
bash collect_and_filter_ssu_pipeline.sh <CONFIG_FILE>
```

Example:

```
bash collect_and_filter_ssu_pipeline.sh config/bacteria_config.conf
```

#### Configuragion file

In configuragion file, input data, working directory and dependencies should be defined by initializing bash environment variables. The following variables are **input-output** veriables:

1. `PREFIX` -- prefix for output files.

2. `WORKDIR` -- working directory, i.e. directory for output and temporary files.

3. `ASSEMBLY_IDS_FPATH` -- file of [Assembly](https://www.ncbi.nlm.nih.gov/assembly/) IDs (one per line).

4. `GENOMES_GBK_DIR` -- directory for `.gbk.gz` files of downloaded genomes.

5. `CONSERVED_REGIONS_FASTA` -- fasta file of the 16S rDNA conserved regions. We used the bundled file `db_creation_and_filtering/conserved_regions.fasta`.

The following variables define **dependencies**:

1. `SEQKIT` -- [seqkit](https://github.com/shenwei356/seqkit) executable.

2. `RANKEDLINEAGE_FPATH` -- file `rankedlineage.dmp` from [NCBI Taxonomy ftp site](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump)

3. `CMSEARCH_FOR_EXTRACT_16S` -- the program `cmsearch` from [Infernal](http://eddylab.org/infernal/). The program is used to extract 16S rRNA genes from genomic sequences. Here, we should use `cmsearcn` 1.1.1, to be fully consistent wuth PGAP.

4. `RFAM_FOR_EXTRACT_16S` -- [Rfam](https://rfam.xfam.org/) covariance model (uncompressed `.cm` file). It is used to extract 16S rRNA genes from genomic sequences. Here, we should use [Rfam 12.0](ftp://ftp.ebi.ac.uk/pub/databases/Rfam/12.0/Rfam.cm.gz), to be fully consistent wuth PGAP.

5. `CMSCAN_FOR_FILTERING` -- the program `cmscan` from [Infernal](http://eddylab.org/infernal/). The program is used to filter incomplete sequences of 16S rRNA genes. Here, we should use the latest version of `cmsearcn` (1.1.4).

6. `CMPRESS_FOR_FILTERING` -- the program `cmpress` from [Infernal](http://eddylab.org/infernal/). The program is used to filter incomplete sequences of 16S rRNA genes. Here, we should use the latest version of `cmpress` (1.1.4).

7. `RFAM_FOR_FILTERING` -- [Rfam](https://rfam.xfam.org/) covariance model (uncompressed `.cm` file). It is used to extract 16S rRNA genes from genomic sequences. Here, we should use [Rfam 14.6](ftp://ftp.ebi.ac.uk/pub/databases/Rfam/14.6/Rfam.cm.gz), to be fully consistent wuth PGAP.

8. `MUSCLE` -- [MUSCLE](https://www.drive5.com/muscle/) executable for multiple sequence alignment.

9. `RFAM_FAMILY_ID` -- ID of specific Rfam covariance model to use. We used the following IDs to predict 16S rRNA: `RF00177` for bacteria and `RF01959` for archaea.


And one variable influences the **behaviour** of the pipeline:

1. `CHECK_CONSERV_REGIONS` -- (accepted values: `0`, `1`). `1` -- check the 16S rDNA conserved regions while finding aberrant genes. `0` -- do not check the conserved regions; detect aberrant genes only by deletions instead.
