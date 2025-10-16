#!/usr/bin/env python3
import sys
from Bio import SeqIO

if len(sys.argv) != 4:
    print("Usage: python grep_gbk.py input.gbk ids.txt output.gbk", file=sys.stderr)
    sys.exit(1)

in_gbk, id_file, out_gbk = sys.argv[1], sys.argv[2], sys.argv[3]

with open(id_file, 'r') as fh:
    wanted = {line.strip() for line in fh if line.strip()}

with open(out_gbk, "w") as gbk_handle:
    for rec in SeqIO.parse(in_gbk, "genbank"):
        if rec.id in wanted:
            SeqIO.write(rec, gbk_handle, "genbank")
