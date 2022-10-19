
import re


def make_seqID(ass_id, acc, seq_start, seq_end, strand_word):
    return 'G_{}:{}:{}-{}:{}'.format(
        ass_id,
        acc,
        seq_start, seq_end,
        strand_word
    )
# end def


def update_seqID(old_seqID, ass_id):
    # From: NZ_CP099994.1:194871-196342_plus
    # To:   G_13242051:NZ_CP100355.1:603895-605369:minus

    old_seqID_pattern = r'([A-Z0-9_\.]+):([0-9]+)-([0-9]+)_([a-z]+)'
    match_obj = re.match(old_seqID_pattern, old_seqID)

    if match_obj is None:
        return old_seqID
    else:
        acc         = match_obj.group(1)
        seq_start   = match_obj.group(2)
        seq_end     = match_obj.group(3)
        strand_word = match_obj.group(4)
        return make_seqID(ass_id, acc, seq_start, seq_end, strand_word)
    # end if
# end def


def update_seqID_column(df, ass_id):

    def set_new_seqID(row):
        old_seqID = row['seqID']
        row['seqID'] = update_seqID(old_seqID, ass_id)
        return row
    # end def

    df = df.apply(set_new_seqID, axis=1)

    return df
# end def
