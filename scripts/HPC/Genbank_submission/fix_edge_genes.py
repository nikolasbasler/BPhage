#!/usr/bin/env python3

# python fix_edge_genes.py genbank_to_table/complete_genomes_fixed.tbl contigs_with_errors_no_stop contigs_with_errors_bad_start genbank_to_table/complete_genomes_fixed_twice.tbl


import sys

in_tbl = sys.argv[1]
missing_stop = sys.argv[2]
bad_start = sys.argv[3]
out_tbl = sys.argv[4]

def is_contig_without_start(id_line):
    contig = id_line.split()[1]
    if contig in no_start:
        return(True)
    else:
        return(False)
    
def is_contig_without_stop(id_line):
    contig = id_line.split()[1]
    if contig in no_stop:
        return(True)
    else:
        return(False)

def contig_length(id_line):
    length = int(id_line.split("_")[1])
    return(length)

with open(missing_stop) as ms, open(bad_start) as bs:
    no_start = [line.rstrip() for line in bs]
    no_stop = [line.rstrip() for line in ms]

with open(in_tbl) as tbl, open(out_tbl, "w") as outfile:
    for line in tbl:
        if line.startswith(">"):
            has_no_start = is_contig_without_start(line)
            has_no_stop = is_contig_without_stop(line)
            cont_length = contig_length(line)
            outfile.write(line)

        elif line.split('\t')[2].strip() == "CDS":
            coords = [line.split('\t')[0], line.split('\t')[1]]

            if has_no_start and coords[0] in ("1","2","3"):
                coords[0] = "<" + coords[0]
            if has_no_start and coords[0] in (str(cont_length), str(cont_length-1), str(cont_length-2)):
                coords[0] = "<" + coords[0]

            if has_no_stop and coords[1] in ("1", "2", "3"):
                coords[1] = ">" + coords[1]
            if has_no_stop and coords[1] in (str(cont_length), str(cont_length-1), str(cont_length-2)):
                coords[1] = ">" + coords[1]

            adapted_line = coords[0] + "\t" + coords[1] + "\tCDS\n"
            outfile.write(adapted_line)
        else:
            outfile.write(line)