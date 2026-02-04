#!/usr/bin/env python3

import os
import sys
import json
import gzip


if len(sys.argv) < 3:
    print('Error!')
    print('Usage:\n  python3{} <FILE.JSON.GZ> <ASSEMBLY_ACCESSION_TO_REMOVE>'.format(sys.argv[0]))
    sys.exit(1)
# end if

json_fpath = sys.argv[1]
asm_acc_to_rm = sys.argv[2]

if not os.path.exists(json_fpath):
    print('Error: file does not exist: `{}`'.format(json_fpath))
    sys.exit(1)
# end if

with gzip.open(json_fpath, 'rt') as input_handle:
    dictionary = json.load(input_handle)
# end with

if not asm_acc_to_rm in dictionary.keys():
    print('Error!')
    print('Cannot find assembly acession `{}` in the json file'.format(asm_acc_to_rm))
    print('Refusing to rewrite')
    sys.exit(0)
# end if

del dictionary[asm_acc_to_rm]

with gzip.open(json_fpath, 'wt') as output_handle:
    json.dump(dictionary, output_handle)
# end with

print('Completed: `{}`'.format(json_fpath))
sys.exit(0)
