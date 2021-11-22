#!/usr/bin/env python3

import os


print(f'\n|=== STARTING SCRIPT `{os.path.basename(__file__)}` ===|\n')

import pandas as pd


# == Input paths ==

# # Bacteria 2.208
# old_ass_acc_fpath = '/mnt/1.5_drive_0/RiboGrove_workdirs/2.208/bacteria/bacteria_refseq_accs_merged.tsv'
# oudated_ass_fpath = '/mnt/1.5_drive_0/RiboGrove_workdirs/2.208/bacteria/outdated_assemblies.tsv'
# ass_status_fpath = '/mnt/1.5_drive_0/RiboGrove_workdirs/2.208/bacteria/assembly_status.tsv'
# outfpath = '/mnt/1.5_drive_0/RiboGrove_workdirs/2.208/bacteria/bacteria_refseq_accs_merged_ass_ids_fixed.tsv'

# Archaea 2.208
old_ass_acc_fpath = '/mnt/1.5_drive_0/RiboGrove_workdirs/2.208/archaea/archaea_refseq_accs_merged.tsv'
oudated_ass_fpath = '/mnt/1.5_drive_0/RiboGrove_workdirs/2.208/archaea/outdated_assemblies.tsv'
ass_status_fpath = '/mnt/1.5_drive_0/RiboGrove_workdirs/2.208/archaea/assembly_status.tsv'
outfpath = '/mnt/1.5_drive_0/RiboGrove_workdirs/2.208/archaea/archaea_refseq_accs_merged_ass_ids_fixed.tsv'

# =================



# Read the input dataframes
old_ass_acc_df = pd.read_csv(old_ass_acc_fpath, sep='\t')
oudated_ass_df = pd.read_csv(oudated_ass_fpath, sep='\t')
ass_status_df = pd.read_csv(ass_status_fpath, sep='\t')


updated_ass_id_df = oudated_ass_df[~pd.isnull(oudated_ass_df['updated_accession'])]

# == Fill the maping dictionary ==
all_ass_ids_in_merged = set(old_ass_acc_df['ass_id'])
all_asmaccs = set(oudated_ass_df['init_accession'])

ass_id_map_dict = dict()

for _, row in updated_ass_id_df.iterrows():
    updated_accession = row['updated_accession']
    if updated_accession in all_asmaccs:
        updated_ass_id = oudated_ass_df.query('init_accession == @updated_accession').iloc[0, 0]
        ass_id_map_dict[updated_ass_id] = row['init_ass_id']
    # end if
# end for

# Found the following updated assemblies:
for k, v in ass_id_map_dict.items():
    print(f'updated: {k}; old: {v}')
# end for


# "Outdate" updated ass_ids

print('Making updated Assembly IDs consistent with the old ones')
print('So that no multiple (old and updated) Assembly IDs will correspond to any genome assembly')

new_ass_acc_df = old_ass_acc_df.copy()

updated_ass_keys_set = set(ass_id_map_dict.keys())

def map_updated_ass_ids_to_old(row):
    
    if row['ass_id'] in updated_ass_keys_set:
        row['ass_id'] = ass_id_map_dict[row['ass_id']]
    # end if
    
    return row
# end def map_updated_ass_ids_to_old

new_ass_acc_df = new_ass_acc_df.apply(map_updated_ass_ids_to_old, axis=1)

print('done\n')


# == Remove non-"Complete genome" and non-"Chromosome" assemblies ==

status_to_rm = {'Contig', 'Scaffold'}
print(f'Removing assemblies of the following completeness levels: {", ".join(status_to_rm)}')

print('Numbers of assemblies having specific completeness level:')
print(ass_status_df.groupby('status', as_index=False).agg({'ass_id': lambda x: x.nunique()}))
print()


# Make everything uppercase
status_to_rm = set(map(str.upper, status_to_rm))
ass_status_df['status'] = ass_status_df['status'].map(str.upper)

ass_ids_to_rm = set(
    ass_status_df.query('status in @status_to_rm')['ass_id']
)

all_asm_count = new_ass_acc_df['ass_id'].nunique()

new_ass_acc_df = new_ass_acc_df.query('not ass_id in @ass_ids_to_rm')

fixed_asm_count = new_ass_acc_df['ass_id'].nunique()

print(f'{all_asm_count-fixed_asm_count} Assembly IDs removed;')
print(f'{fixed_asm_count} Assembly IDs kept;')

new_ass_acc_df.to_csv(
    outfpath,
    sep='\t',
    header=True,
    index=False,
    na_rep='NA'
)

print(outfpath)
print('Completed!')

print(f'\n|=== EXITTING SCRIPT `{os.path.basename(__file__)}` ===|\n')

