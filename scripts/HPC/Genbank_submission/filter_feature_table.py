#!/usr/bin/env python3
import sys

in_tbl  = sys.argv[1]
id_file = sys.argv[2]
out_tbl = sys.argv[3]

wanted = {line.strip() for line in open(id_file) if line.strip()}

with open(in_tbl, "r") as tbl, open(out_tbl, "w") as outfile:
    is_in_wanted = False
    for line in tbl:
        if line.startswith(">Feature"):
            parts = line.strip().split()
            # guard against malformed lines
            if len(parts) < 2:
                is_in_wanted = False
                continue
            contig = parts[1]
            if contig in wanted:
                is_in_wanted = True
                outfile.write(line)     # write the header for a wanted block
            else:
                is_in_wanted = False    # stop writing at this header (don't write it)
        elif is_in_wanted:
            outfile.write(line)        # write non-header lines only when we're inside a wanted block
