# RiboGrove creation demo

Here is the example of a RiboGrove working directory and of a `.conf` file alogn with file `demo_assembly_ids.txt` listing three Assembly Ids.

To test this demo, you first should met all dependencies for the pipeline `collect_and_filter_ssu_pipeline.sh`. Then, you should modify paths in file `demo.conf`: all input files and files of dependencies should exists, and output directories (`WORKDIR` and `GENOMES_GBK_DIR`) should satisfy you.

Then, you can run the following command:

```bash
cd demo/
bash ../db_creation_and_filtering/collect_and_filter_ssu_pipeline.sh demo.conf
```
