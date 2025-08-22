#!/usr/bin/env python3
import sys
from Bio import SeqIO

in_gbk  = sys.argv[1]
id_file = sys.argv[2]
out_gbk = sys.argv[3]

wanted = {line.strip() for line in open(id_file) if line.strip()}

with open(out_gbk, "w") as gbk_handle:
    for rec in SeqIO.parse(in_gbk, "genbank"):
        if rec.id in wanted:
            SeqIO.write(rec, gbk_handle, "genbank")
