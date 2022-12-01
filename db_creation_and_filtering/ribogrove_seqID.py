
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
    # It seqiD is new, we'll just update ass_id:
    # From: G_7144461:NZ_CP053223.1:191894-193388:plus
    # To:   G_14020671:NZ_CP053223.1:191894-193388:plus
    #         ~~~~~~~~

    old_seqID_pattern = r'([A-Z0-9_\.]+):([0-9]+)-([0-9]+)_([a-z]+)'
    old_match = re.match(old_seqID_pattern, old_seqID)

    if not old_match is None:
        acc         = old_match.group(1)
        seq_start   = old_match.group(2)
        seq_end     = old_match.group(3)
        strand_word = old_match.group(4)
        return make_seqID(ass_id, acc, seq_start, seq_end, strand_word)
    else:
        new_seqID_pattern = r'G_([0-9]+):([A-Z0-9_\.]+):([0-9]+)-([0-9]+):([a-z]+)'
        new_match = re.match(new_seqID_pattern, old_seqID)
        if not new_match is None:
            old_ass_id = new_match.group(1)
            new_seq_id = old_seqID.replace(
                'G_{}:'.format(old_ass_id),
                'G_{}:'.format(ass_id),
            )
            return new_seq_id
        else:
            err_msg = '\nError: cannot recognize seqID: `{}`'.format(old_seqID)
            raise ValueError(err_msg)
        # end if
        return old_seqID
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
