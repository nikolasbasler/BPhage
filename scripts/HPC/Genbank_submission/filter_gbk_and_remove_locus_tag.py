#!/usr/bin/env python3
import sys
from Bio import SeqIO

in_gbk  = sys.argv[1]
id_file = sys.argv[2]
out_gbk = sys.argv[3]

with open(id_file) as fh:
    wanted = {line.strip() for line in fh if line.strip()}

with open(out_gbk, "w") as gbk_handle:
    for rec in SeqIO.parse(in_gbk, "genbank"):
        if rec.id in wanted:
            # Remove any locus_tag qualifiers from every feature (case-insensitive)
            for feat in rec.features:
                # list() to avoid runtime-dict-change issues
                remove_keys = [k for k in list(feat.qualifiers.keys()) if k.lower() == "locus_tag"]
                for k in remove_keys:
                    feat.qualifiers.pop(k, None)
            SeqIO.write(rec, gbk_handle, "genbank")

