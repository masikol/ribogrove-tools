
# Create RiboGrove

The process of RiboGrove release creation consists of two major steps:

1. `collect_and_filter/` -- collect 16S rRNA gene sequences and filter them. This must be done separately for bacteria and archaea.

2. `combine_release/` -- combine bacterial and archaeal data into a RiboGrove release, namely into files `ribogrove_XX.XXX_sequences.fasta.gz` and `metadata.zip` (`XX.XXX` is RiboGrove release number, e.g. 10.216).
