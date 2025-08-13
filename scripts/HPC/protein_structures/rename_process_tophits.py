#!/usr/bin/env python3
import argparse
from argparse import RawTextHelpFormatter
from Bio import SeqIO
import os
import shutil
import sys
import json
from pathlib import Path
import numpy as np
from Bio.Seq import Seq
import pandas as pd


def get_input():
    usage = "python3 rename_preocess_tophits.py ..."
    parser = argparse.ArgumentParser(
        description="script to rename top ranked pdbs and json from colabfold output (e.g. 0_unrelaxed_rank_001_alphafold2_ptm_model_1_seed_000.pdb) to the name as per its multifasta and get some more metadata.",
        formatter_class=RawTextHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--infile",
        action="store",
        help="Input file of AAs. Headers with the names",
        required=True,
    )
    parser.add_argument(
        "-p", "--predictions_dir", action="store", help="PDB Dir", required=True
    )
    parser.add_argument(
        "-o", "--outdir", action="store", help="Output directory", required=True
    )
    parser.add_argument(
        "-f", "--force", help="Overwrites the output directory.", action="store_true"
    )
    args = parser.parse_args()

    return args


def main():
    args = get_input()

    if args.force == True:
        if os.path.isdir(args.outdir) == True:
            print(
                f"Removing output directory {args.outdir} as -f or --force was specified."
            )
            shutil.rmtree(args.outdir)
        elif os.path.isfile(args.outdir) == True:
            os.remove(args.outdir)
        else:
            print(
                f"--force was specified even though the output directory {args.outdir} does not already exist. Continuing."
            )
    else:
        if os.path.isdir(args.outdir) == True or os.path.isfile(args.outdir) == True:
            print(
                f"The output directory {args.outdir} already exists and force was not specified. Please specify -f or --force to overwrite it."
            )
            sys.exit()

    # make outdir if it doesn't exist
    os.makedirs(args.outdir, exist_ok=True)
    pdb_out_dir = os.path.join(args.outdir, "renamed_pdbs")
    json_out_dir = os.path.join(args.outdir, "renamed_jsons")
    pae_out_dir = os.path.join(args.outdir, "renamed_pae_jsons")

    os.makedirs(pdb_out_dir, exist_ok=True)
    os.makedirs(json_out_dir, exist_ok=True)
    os.makedirs(pae_out_dir, exist_ok=True)

    # Initialize an empty list to store protein headers
    protein_headers = []

    protein_len_dict = {}

    for record in SeqIO.parse(args.infile, "fasta"):
        header = f"{record.id}"
        protein_headers.append(header)
        protein_len_dict[header] = {"length": len(Seq(record.seq))}

    # List all files in the directory ending with ".pdb"
    pdb_files = [f for f in os.listdir(args.predictions_dir) if f.endswith(".pdb") and "rank_001" in f]

    plddt_ptm_dict = {}

    # Iterate through each PDB file
    for pdb_file in pdb_files:
        # Split the filename by "_" to get the index
        file_parts = pdb_file.split("_")

        # Extract the index from the first part of the filename
        index = int(file_parts[0])

        # Get the corresponding protein header from the list
        if 0 <= index < len(protein_headers):
            new_name = protein_headers[index] + ".pdb"

            new_path = os.path.join(pdb_out_dir, new_name)
            shutil.copy(os.path.join(args.predictions_dir, pdb_file), new_path)
            
            # json

            json_input_files = [
                os.path.join(
                    args.predictions_dir,
                    f"{index}_scores_rank_001_alphafold2_ptm_model_1_seed_000.json",
                ),
                os.path.join(
                    args.predictions_dir,
                    f"{index}_scores_rank_001_alphafold2_ptm_model_2_seed_000.json",
                ),
                os.path.join(
                    args.predictions_dir,
                    f"{index}_scores_rank_001_alphafold2_ptm_model_3_seed_000.json",
                ),
            ]

            # check if any exist
            if (
                os.path.exists(json_input_files[0])
                or os.path.exists(json_input_files[1])
                or os.path.exists(json_input_files[2])
            ):
                new_name_json = protein_headers[index] + ".json"
                new_path_json = os.path.join(json_out_dir, new_name_json)
                if os.path.exists(json_input_files[0]):
                    input_json = json_input_files[0]
                elif os.path.exists(json_input_files[1]):
                    input_json = json_input_files[1]
                elif os.path.exists(json_input_files[2]):
                    input_json = json_input_files[2]
                # copy the json
                shutil.copy(input_json, new_path_json)
                # get the data
                # parse the json
                scores = json.loads(Path(new_path_json).read_text())

                # get the stats
                plddt = np.asarray(scores["plddt"])
                mean_plddt = round(np.mean(plddt), 2)
                ptm = scores["ptm"]

                # Add values to the dictionary
                plddt_ptm_dict[protein_headers[index]] = {
                    "folded": True,
                    "mean_plddt": mean_plddt,
                    "ptm": ptm,
                }

                # 
                pae_file = f"{index}_predicted_aligned_error_v1.json"
                new_name_pae = protein_headers[index] + "_pae.json"
                new_path_pae = os.path.join(os.path.join(pae_out_dir, new_name_pae))
                shutil.copy(os.path.join(args.predictions_dir, pae_file), new_path_pae)

            else:
                print(
                    f"{file_parts} aka {new_name} json does not exist. Please copy it in"
                )

                folded = False

                # Add values to the dictionary
                plddt_ptm_dict[protein_headers[index]] = {
                    "folded": False,
                    "mean_plddt": mean_plddt,
                    "ptm": ptm,
                }

        else:
            print(f"Skipping {pdb_file} as index is out of range.")

        # Convert dictionary to pandas DataFrame
    len_df = pd.DataFrame.from_dict(protein_len_dict, orient="index")
    len_df = len_df.rename_axis("protein")

    plddt_df = pd.DataFrame.from_dict(plddt_ptm_dict, orient="index")
    plddt_df = plddt_df.rename_axis("protein")

    # Merge the DataFrames on the 'protein' column
    merged_df = pd.merge(len_df, plddt_df, on="protein", how="inner")

    # Save DataFrame to a CSV file
    outfile_path = os.path.join(args.outdir, "plddt_ptm_len.tsv")
    merged_df.to_csv(outfile_path, sep="\t", index=True)


if __name__ == "__main__":
    main()
