#!/usr/bin/env python3
import sys
from Bio import SeqIO

in_gbk  = sys.argv[1]
id_file = sys.argv[2]

wanted = {line.strip() for line in open(id_file) if line.strip()}

with open(in_gbk, "r") as gbk_handle:
    for rec in SeqIO.parse(gbk_handle, "genbank"):
        if rec.id in wanted:
            print(rec.id)
            # Print out first and last CDS entry
            cds_features = [feature for feature in rec.features if feature.type == "CDS"]
            if cds_features:
                print("First CDS:", cds_features[0])
                print("Last CDS:", cds_features[-1])