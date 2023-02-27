
# SeqID versions:
# v1: NZ_CP099994.1:194871-196342_plus
# v2: G_13242051:NZ_CP100355.1:603895-605369:minus
# v3: GCF_024204665.1:NZ_CP100355.1:603895-605369:minus
# v3 is the latest and the current

import re


V1_PATTERN = re.compile(
    r'([A-Z0-9_\.]+):([0-9]+)-([0-9]+)_([a-z]+)'
)
V2_PATTERN = re.compile(
    r'G_([0-9]+):([A-Z0-9_\.]+):([0-9]+)-([0-9]+):([a-z]+)'
)
V3_PATTERN = re.compile(
    r'(GCF_[0-9\.]+):([A-Z0-9_\.]+):([0-9]+)-([0-9]+):([a-z]+)'
)


def make_seqID(asm_acc, acc, seq_start, seq_end, strand_word):
    return '{}:{}:{}-{}:{}'.format(
        asm_acc,
        acc,
        seq_start, seq_end,
        strand_word
    )
# end def

def make_seqID_v2(ass_id, acc, seq_start, seq_end, strand_word):
    return 'G_{}:{}:{}-{}:{}'.format(
        ass_id,
        acc,
        seq_start, seq_end,
        strand_word
    )
# end def


def update_seqID_v2_to_v3(old_seqID, asm_acc):
    # From: G_13242051:NZ_CP100355.1:603895-605369:minus
    # To:   GCF_024204665.1:NZ_CP100355.1:603895-605369:minus

    global V2_PATTERN
    global V3_PATTERN

    old_match = re.match(V2_PATTERN, old_seqID)

    if not old_match is None:
        seq_acc     = old_match.group(2)
        seq_start   = old_match.group(3)
        seq_end     = old_match.group(4)
        strand_word = old_match.group(5)
        return make_seqID(asm_acc, seq_acc, seq_start, seq_end, strand_word)
    else:
        new_match = re.match(V3_PATTERN, old_seqID)
        if new_match is None:
            err_msg = '\nError: cannot recognize seqID: `{}`'.format(old_seqID)
            raise ValueError(err_msg)
        # end if
        return old_seqID
    # end if
# end def


def update_seqID_column_v2_to_v3(df, asm_acc):

    def set_new_seqID(row):
        old_seqID = row['seqID']
        row['seqID'] = update_seqID_v2_to_v3(old_seqID, asm_acc)
        return row
    # end def

    df = df.apply(set_new_seqID, axis=1)

    return df
# end def


def parse_asm_acc(seqID):
    return seqID.partition(':')[0]
# end def

def parse_ass_id_v2(seqID):
    return seqID.partition(':')[0][2:]
# end def

def parse_seq_acc(seqID):
    return seqID.split(':')[1]
# end def



def update_seqID_v1_to_v2(old_seqID, ass_id):
    # From: NZ_CP099994.1:194871-196342_plus
    # To:   G_13242051:NZ_CP100355.1:603895-605369:minus
    # It seqiD is new, we'll just update ass_id:
    # From: G_7144461:NZ_CP053223.1:191894-193388:plus
    # To:   G_14020671:NZ_CP053223.1:191894-193388:plus
    #         ~~~~~~~~
    global V1_PATTERN
    global V2_PATTERN

    old_match = re.match(V1_PATTERN, old_seqID)

    if not old_match is None:
        acc         = old_match.group(1)
        seq_start   = old_match.group(2)
        seq_end     = old_match.group(3)
        strand_word = old_match.group(4)
        return make_seqID(ass_id, acc, seq_start, seq_end, strand_word)
    else:
        new_match = re.match(V2_PATTERN, old_seqID)
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


def update_seqID_column_v1_to_v2(df, ass_id):

    def set_new_seqID(row):
        old_seqID = row['seqID']
        row['seqID'] = update_seqID_v1_to_v2(old_seqID, ass_id)
        return row
    # end def

    df = df.apply(set_new_seqID, axis=1)

    return df
# end def
