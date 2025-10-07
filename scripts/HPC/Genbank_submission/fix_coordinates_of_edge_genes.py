#!/usr/bin/env python3
# python $repo_location/scripts/HPC/Genbank_submission/fix_coordinates_of_edge_genes.py temp.fixed_repeat_regions temp.fixed_clipped_cds

import sys

in_tbl = sys.argv[1]
out_tbl = sys.argv[2]

with open(in_tbl) as tbl, open(out_tbl, "w") as outfile:
    for line in tbl:
        if line.startswith(">"):
            cont_length = line.split("_")[1]
            outfile.write(line)

            id = line.split()[1].strip()

        elif line.split('\t')[2].strip() == "CDS":
            coords = [line.split('\t')[0], line.split('\t')[1]]
            stripped_coords = [int(c.lstrip("<>").rstrip(">")) for c in coords]

            if coords == [str(stripped_coords[0]), str(stripped_coords[1])]:
                frameline = "\t\t\tcodon_start\t1\n"
                outfile.write(line)
                outfile.write(frameline)
                continue

            else:
                three_prime_is_clipped = coords[1].startswith(">")
                five_prime_is_clipped = coords[0].startswith("<")

                if three_prime_is_clipped and stripped_coords[0] < stripped_coords[1]:
                    # Forward CDS, 3' clipped
                    outline = coords[0] + "\t>" + cont_length + "\tCDS\n"
                    frameline = "\t\t\tcodon_start\t1\n"
                    outfile.write(outline)
                    outfile.write(frameline)
                    continue
                if three_prime_is_clipped and stripped_coords[0] > stripped_coords[1]:
                    # Reverse CDS, 3' clipped
                    outline = coords[0] + "\t>" + "1" + "\tCDS\n"
                    frameline = "\t\t\tcodon_start\t1\n"
                    outfile.write(outline)
                    outfile.write(frameline)
                    continue
                # if five_prime_is_clipped and stripped_coords[0] < stripped_coords[1] and not stripped_coords[0] == 1:
                if five_prime_is_clipped and stripped_coords[0] < stripped_coords[1]:
                    # Forward CDS, 5' clipped
                    outline = "<1" + "\t" + coords[1] + "\tCDS\n"
                    outfile.write(outline)
                    frameline = "\t\t\tcodon_start\t" + str(stripped_coords[0]) + "\n"
                    outfile.write(frameline)
                    continue
                # if five_prime_is_clipped and stripped_coords[0] > stripped_coords[1] and not stripped_coords[0] == int(cont_length):
                if five_prime_is_clipped and stripped_coords[0] > stripped_coords[1]:
                    # Reverse CDS, 5' clipped
                    outline = "<" + cont_length + "\t" + coords[1] + "\tCDS\n"
                    outfile.write(outline)

                    frame = str(int(cont_length) - int(stripped_coords[0]) + 1)
                    frameline = "\t\t\tcodon_start\t" + frame + "\n"
                    outfile.write(frameline)
                    
                    continue
                else:
                    outfile.write(outline)
        else:
            outfile.write(line)
