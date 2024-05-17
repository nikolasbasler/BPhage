#!/usr/bin/env python3

# 1st argument: Input fasta file
# 2nd argument: Max proportion of Xs allowed per sequence. If set to 0, proportion of Xs per sequence is printed instead.
# 
# Example:
# python3 remove_seqs_with_Xs_from_fasta input.fasta 0.5

import sys
from Bio import SeqIO

if float(sys.argv[2])==0:
    for record in SeqIO.parse(sys.argv[1], "fasta"):
        print(record.id + "\t" + str(record.seq.count('X') / len(record.seq)))

else:
    for record in SeqIO.parse(sys.argv[1], "fasta"):
        if record.seq.count('X') / len(record.seq) <= float(sys.argv[2]):
            print(record.format("fasta"), end='')
