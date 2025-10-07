#!/usr/bin/env python3

# python $repo_location/scripts/HPC/Genbank_submission/fix_tRNAs.py temp.fixed_clipped_cds trna_anticodons_and_aminoacids temp.fixed_trnas

import sys

in_tbl = sys.argv[1]
trna_anticodons = sys.argv[2]
out_tbl = sys.argv[3]

# def rev_comp_RNA(codon):
#     # This doesn't take wobble into account... So I'm not using if for now. But according to the NLM support, it's not needed anyway...
#     complement = {"A": "U", "U": "A", "C": "G", "G": "C", "N": "N"}
#     return "".join(complement[base] for base in reversed(codon))

contigs_with_trna=[]
trna_loci = {}
with open(trna_anticodons) as ta:
    for line in ta:
        contig = line.split()[0].strip()
        coord_one = line.split()[1].strip()
        coord_two = line.split()[2].strip()
        anticodon = line.split()[3].strip()
        # rev_anticodon = rev_comp_RNA(anticodon)
        aminoacid = line.split()[4].strip()

        contigs_with_trna.append(contig)

        locus = "+".join([contig, coord_one, coord_two])
        # trna_loci[locus] = (rev_anticodon, aminoacid)
        trna_loci[locus] = (anticodon, aminoacid)

addition_lines = 0
with open(in_tbl) as tbl, open(out_tbl, "w") as outfile:
    for line in tbl:
        if line.startswith(">"):
            contig_name=line.split()[1].strip()
            outfile.write(line)
            continue
        if contig_name in contigs_with_trna and line.split('\t')[2].strip() == "tRNA":
            locus = "+".join([contig_name, line.split('\t')[0].strip(), line.split('\t')[1].strip()])
            anticodon = trna_loci[locus][0]
            aminoacid = trna_loci[locus][1]

            noteline = "\t\t\tnote\tisotype=" + aminoacid + ":anticodon=" + anticodon
            if aminoacid == "Undet":
                noteline = "\t\t\tnote\tundetermined amino acid"
                aminoacid = "OTHER"
            if aminoacid == "Sup":
                noteline = noteline + ":suppressor tRNA"
                aminoacid = "OTHER"

            product_outline = "\t\t\tproduct\ttRNA-" + aminoacid + "\n"

            outfile.write(line)
            outfile.write(product_outline)
            outfile.write(noteline + "\n")

            addition_lines += 2
        else:
            outfile.write(line)
print(addition_lines, "additional lines were added.")
