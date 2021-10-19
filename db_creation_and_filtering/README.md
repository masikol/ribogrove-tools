# RybaSom database creation and sequence filtering

Here you can find scripts that were (and can be) used for RybaSom construction.

Each Python script can be used on it's own using command line interace. At the same time, they are integrated in a pipeline `collect_and_filter_ssu_pipeline.sh`.

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

In configuragion file, input data, working directory and dependencies should be defined by initializing bash environment variables. Necessary variables are following:

1. `PREFIX`: prefix for output files. 

1. File of [Assembly](https://www.ncbi.nlm.nih.gov/assembly/) IDs (one per line).

2. 


