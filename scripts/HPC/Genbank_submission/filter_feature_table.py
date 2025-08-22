#!/usr/bin/env python3
import sys

in_tbl  = sys.argv[1]
out_tbl = sys.argv[2]

wanted = {line.strip() for line in open(id_file) if line.strip()}

with open(in_tbl, "r") as tbl, open(out_tbl, "w") as outfile:
    is_in_wanted = False
    for line in tbl:
        if line.startswith(">Feature") and line.strip().split()[1] in wanted:
            is_in_wanted = True
            outfile.write("contig_" + line.split()[1] + "\n")
            continue
        if is_in_wanted:
            outfile.write(line)
        if line.startswith(">Feature") and line.strip().split()[1] not in wanted:
            is_in_wanted = False
