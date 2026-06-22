
import hashlib


def make_asm_seqID_hash(seqIDs: list[str]) -> str:
    sep = '$'
    seqIDs_joint = sep.join(sorted(seqIDs))
    return hashlib.md5(seqIDs_joint.encode()).hexdigest()
# end def
