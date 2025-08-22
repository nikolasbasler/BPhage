#!/usr/bin/env python3
import sys

in_tbl  = sys.argv[1]
id_file = sys.argv[2]
out_tbl = sys.argv[3]
out_crcl = sys.argv[4]

wanted = {line.strip() for line in open(id_file) if line.strip()}

def is_cutoff_at_beginning(coord_line):
    if coord_line.split('\t')[0].startswith("<"):
        return(True)
    else:
        return(False)

def is_cutoff_at_end(coord_line):
    if coord_line.split('\t')[1].startswith(">"):
        return(True)
    else:
        return(False)
    
def get_coords(coord_line):
    if is_cutoff_at_beginning(coord_line):
        start = int(coord_line.split('\t')[0][1:])
    else:
        start = int(coord_line.split('\t')[0])
    if is_cutoff_at_end(coord_line):
        end = int(coord_line.split('\t')[1][1:])
    else:
        end = int(coord_line.split('\t')[1])
    return (start, end)

def is_reverse(coord_line):
    coords = get_coords(coord_line)
    if coords[0] > coords[1]:
        return(True)
    else:
        return(False)

def length_divisible_by_three(coord_line):
    coords = get_coords(coord_line)
    length_of_seq = abs(coords[1] - coords[0]) + 1
    if length_of_seq % 3 == 0:
        return(True)
    else:
        return(False)
    
def get_contig_length(contig):
    return int(contig.split("_")[1])

with open(in_tbl, "r") as tbl, open(out_tbl, "w") as outfile:
# with open(in_tbl, "r") as tbl:
    is_in_wanted = False
    new_contig = False
    pseudo_circle = {}
    for line in tbl:
        if line.startswith(">Feature") and line.strip().split()[1] in wanted:
            outfile.write("contig_" + line.split()[1] + "\n")
            # print(line.strip().split()[1])
            is_in_wanted = True
            contig = line.strip().split()[1]
            new_contig = True
            continue
        if line.startswith(">Feature") and line.strip().split()[1] not in wanted:
            is_in_wanted = False
            continue
        if is_in_wanted and not line.startswith("\t"):
            outline = line.strip() + \
                " " + contig + \
                " cutoff_at_beginning:" + str(is_cutoff_at_beginning(line.strip())) + \
                " cutoff_at_end:" + str(is_cutoff_at_end(line.strip())) + \
                " reverse:" + str(is_reverse(line.strip())) + \
                " length_divisible_by_three:" + str(length_divisible_by_three(line.strip()))
            # print(outline)
            outfile.write(outline + "\n")
            # if new_contig:
            #     pseudo_circle[contig + "_first"] = get_coords(line.strip())
            #     new_contig = False
            # if not new_contig:
            #     pseudo_circle[contig + "_last"] = get_coords(line.strip())

            if new_contig and is_cutoff_at_beginning(line.strip()):
                pseudo_circle[contig] = get_contig_length(contig), get_coords(line.strip()), (None, None)
                new_contig = False
            if not new_contig and is_cutoff_at_end(line.strip()):
                pseudo_circle[contig] = get_contig_length(contig), pseudo_circle[contig][1], get_coords(line.strip())


with open(out_crcl, "w") as crcl:
    for contig in pseudo_circle:
        if pseudo_circle[contig][2] == (None, None):
            continue
        if pseudo_circle[contig][1][0] > pseudo_circle[contig][1][1]:
            first_cds_is_reverse = True
        else:
            first_cds_is_reverse = False
        if pseudo_circle[contig][2][0] > pseudo_circle[contig][2][1]:
            last_cds_is_reverse = True
        else:
            last_cds_is_reverse = False

        if first_cds_is_reverse != last_cds_is_reverse:
            # print(contig + ": First and last CDS not in same direction.")
            crcl.write(contig + ": First and last CDS not in same direction.\n")
            continue

        contig_length = get_contig_length(contig)

        length_of_first_cds = max(pseudo_circle[contig][1])
        length_of_last_cds = contig_length - min(pseudo_circle[contig][2]) + 1

        if (length_of_first_cds + length_of_last_cds) % 3 != 0:
            # print(contig + ": First and last CDS lengths are NOT divisible by 3.")
            crcl.write(contig + ": Not divisible by 3 if circularized.\n")
