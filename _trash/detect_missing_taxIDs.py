#!/usr/bin/env python3

rankedlineage_path = '/mnt/1.5_drive_0/16S_scrubbling/taxonomy/new_taxdump/ranked_lineage_just_tabs.dmp'
my_taxIDs_fpath = '/mnt/1.5_drive_0/16S_scrubbling/taxonomy/taxIDs.tsv'

all_seqIDs = {int(line.split('\t')[0]) for line in open(rankedlineage_path).readlines()}

print('Missing taxIDs:')

with open(my_taxIDs_fpath, 'rt') as my_taxIDs_fpath:
    print(my_taxIDs_fpath.readline().strip()) # pass header
    for line in my_taxIDs_fpath:
        taxID = int(line.split('\t')[2])
        if taxID not in all_seqIDs:
            print(line.strip())
        # end if
    # end for
# wnd with


print('Completed!')
