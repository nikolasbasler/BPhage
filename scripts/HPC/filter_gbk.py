import sys
from Bio import SeqIO

# Will filter out genomes ("LOCUS" entries) from a .gbk file and writes the
# filtered .gbk file to stdout.

# Usage: python filter_gbk.py <input_file.gbk> <identifier1> [<identifier2> ...]

input_file = sys.argv[1]
identifiers_to_remove = set(sys.argv[2:])

for record in SeqIO.parse(input_file, "genbank"):
    if record.id not in identifiers_to_remove:
        SeqIO.write(record, sys.stdout, "genbank")