# RiboGrove creation demo

Here is the example of a RiboGrove working directory and of a `.conf` file.

To test this demo, you first should met all dependencies for the pipeline `collect_and_filter.sh`. Then, you should modify paths in file `demo.conf`: all input files and files of dependencies should exists.

Then, you can run the pipeline in test mode (set `-t` option):

```bash
cd create_RiboGrove/collect_and_filter/

bash collect_and_filter.sh demo/demo.conf -t
```
