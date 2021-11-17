#!/usr/bin/env python3

import sys

subject_infpath = '/mnt/1.5_drive_0/16S_scrubbling/ribogrove-paper/db_creation_and_filtering/accessions_after_title_map.txt'
query_infpath = '/mnt/1.5_drive_0/RiboGrove_workdirs/2.208/bacteria/accs_208_bacteria_noWGS.txt'

subject_set = set(
    map(
        str.strip,
        open(subject_infpath).readlines()
    )
)

with open(query_infpath, 'rt') as q_infile:
    for acc in q_infile:
        acc = acc.strip()
        if not acc in subject_set:
            print(acc)
        # end if
    # end for
# end with

sys.exit(0)
