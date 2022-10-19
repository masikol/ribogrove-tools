
import pandas as pd

from ribogrove_seqID import update_seqID


def update_tblout_file(tblout_fpath, categories_fpath):

    tblout_df = pd.read_csv(tblout_fpath, sep='\t')
    categories_df = pd.read_csv(categories_fpath, sep='\t')

    seqID_assID_dict = {
        row['seqID']: row['ass_id'] for i, row in categories_df.iterrows()
    }
    del categories_df

    def set_new_query_name(row):
        old_seqID = row['query_name']
        ass_id = seqID_assID_dict[old_seqID]
        row['query_name'] = update_seqID(old_seqID, ass_id)
        return row
    # end def

    new_tblout_df = tblout_df.apply(set_new_query_name, axis=1)

    return new_tblout_df
# end def
