import sys
from Bio import SeqIO

# Will filter out genomes ("LOCUS" entries) from a .gbk file and writes the
# filtered .gbk file to stdout.

# Usage: python filter_gbk.py <input_file.gbk> <identifier1> [<identifier2> ...]

input_file = sys.argv[1]
patterns_file = sys.argv[2]
identifiers_to_remove = set(sys.argv[3:])

with open(patterns_file, 'r') as pf:
    patterns_to_remove = [line.strip() for line in pf]

for record in SeqIO.parse(input_file, "genbank"):
    if record.id not in identifiers_to_remove and record.id not in patterns_to_remove:
        SeqIO.write(record, sys.stdout, "genbank")


# input_file = sys.argv[1]
# identifiers_to_remove = set(sys.argv[2:])

# with open(patterns_file, 'r') as pf:
#     patterns_to_remove = [line.strip() for line in pf]

# for record in SeqIO.parse(input_file, "genbank"):
#     if record.id not in identifiers_to_remove:
#         SeqIO.write(record, sys.stdout, "genbank")