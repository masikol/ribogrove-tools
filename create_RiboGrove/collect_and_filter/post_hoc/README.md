This directory contains scripts which can remove all traces a genome assembly from a RiboGrove workdir (e.g. `RiboGrove_23.229/bacteria`).

**Attention!** This scripts do not add assembly accessions to the blacklist. Please, do it yourself.

The main script is `post_hoc_remove_genome.sh`. Usage:

```bash
bash post_hoc_remove_genome.sh \
    <CONFIG_FILE> \
    <ASSEMBLY_ACCESSION_TO_REMOVE>
```

Where `<CONFIG_FILE>` is the same configuration file which must be provided to the script `collect_and_filter.sh`.

Example:

```bash
bash post_hoc_remove_genome.sh \
    RiboGrove_workdirs/23.229/bacteria/bacteria.conf \
    GCF_047323925.1
```
