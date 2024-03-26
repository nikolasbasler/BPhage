#!/usr/bin/env python3
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr

# Argument parser
parser = argparse.ArgumentParser(description="Correlation analysis script.")
parser.add_argument(
    "-i",
    "--input",
    type=str,
    metavar="FILE",
    required=True,
    help="Abundance table file (tsv).",
)
parser.add_argument(
    "-o",
    "--output",
    type=str,
    metavar="OUTPUT",
    required=True,
    help="Prefix for the output name.",
)
parser.add_argument(
    "-s",
    "--segments",
    type=str,
    metavar="FILE",
    help="File with a list of contigs of interest (often RdRP segments), each on a new line.",
)
parser.add_argument(
    "-l",
    "--lengths",
    type=str,
    metavar="FILE",
    help="File with the lengths of each contig.",
)
parser.add_argument(
    "-p",
    "--prevalence",
    type=float,
    metavar="FLOAT",
    default=0.1,
    help="Minimum percentage of samples for correlation analysis. (default: %(default)s)",
)
parser.add_argument(
    "-c",
    "--correlation",
    type=float,
    metavar="FLOAT",
    default=0.5,
    help="Minimum correlation to keep pairs. (default: %(default)s)",
)
parser.add_argument(
    "--strict",
    dest="strict",
    action="store_true",
    default=False,
    help="The correlation threshold should be met for all provided segments. (default: %(default)s)",
)
parser.add_argument(
    "-v",
    "--p_value",
    type=float,
    metavar="FLOAT",
    default=0.05,
    help="Significance threshold for correlations. (default: %(default)s)",
)

def clean_abundance(df):
    df_clean = df.drop(df.columns[df.columns.str.contains("NC")], axis=1)
    df_clean["sample_count"] = df_clean.apply(lambda row: row[row != 0].count(), axis=1)
    df_clean["proportion_samples"] = df_clean["sample_count"] / (df_clean.shape[1] - 1)

    return df_clean

def spearmanr_pval(x,y):
    return spearmanr(x,y)[1]

def pearsonr_pval(x,y):
    return pearsonr(x,y)[1]

def create_correlation_matrix(df_transposed, p_value):
    correlation_matrix = df_transposed.corr(method="spearman")
    p_value_matrix = df_transposed.corr(method=spearmanr_pval)
    # correlation_matrix = df_transposed.corr(method="pearson")
    # p_value_matrix = df_transposed.corr(method=pearsonr_pval)
    correlation_matrix[p_value_matrix > p_value] = np.nan

    mask = np.triu(np.ones(correlation_matrix.shape), k=1).astype(bool)
    masked_correlation_matrix = correlation_matrix.mask(mask)

    correlation_matrix = masked_correlation_matrix.rename_axis(
        axis=0, mapper="Contig1"
    ).rename_axis(axis=1, mapper="Contig2")

    return correlation_matrix


def segment_correlation_matrix(df, segment_list):
    correlation_results_df = pd.DataFrame()

    # Loop through each segment in segment_list
    for i in segment_list:
        df2 = df.loc[i]
        df3 = df.corrwith(df2, axis=1, method="spearman")
        # To do: Implement p-value threshold

        # Append the results to the DataFrame with i as the column name
        correlation_results_df[i] = df3

    return correlation_results_df


def create_segment_list(segment_file):
    file_path = Path(segment_file)

    # Read the file into a list
    with open(file_path, "r") as file:
        segment_list = file.readlines()

    segment_list = [line.strip() for line in segment_list]
    return segment_list


def main():
    args = parser.parse_args()

    print("Read in abundance table.")
    abundance_df = pd.read_csv(args.input, sep="\t", index_col=0)
    df = clean_abundance(abundance_df)

    # Define the threshold
    prevalence_threshold = args.prevalence

    # Filter rows where the proportion of 0s is less than or equal to the threshold
    filtered_df = df[df["proportion_samples"] >= prevalence_threshold].drop(
        ["sample_count", "proportion_samples"], axis=1
    )

    n = len(filtered_df)

    if args.lengths:
        print(f"Using contig length corrected abundance table.")
        lengths = pd.read_csv(
            args.lengths, sep="\t", index_col=0, header=None, names=["Contig", "length"]
        )
        df = filtered_df.div(lengths["length"], axis=0).dropna(how="all")
    else:
        print(f"Using absence/presence abundance table.")
        df = filtered_df.map(lambda x: 1 if x > 0 else x)

    print(
        f"Calculate correlation matrix for {n} contigs (contig prevalence in samples set to {prevalence_threshold*100}%)."
    )

    cor_threshold = args.correlation

    if args.segments:
        segment_list = create_segment_list(args.segments)

        # Check if all values in segment_list are present in df indices
        missing_segments = [
            segment for segment in segment_list if segment not in df.index
        ]

        if missing_segments:
            missing_segments_str = ", ".join(map(str, missing_segments))
            raise ValueError(
                f"The following segment(s) are not present in the sample prevalence filtered abundance table: {missing_segments_str}\n"
                f"Consider lowering the prevalence threshold (-p) which is currently at {prevalence_threshold}"
            )

        correlation_results_df = segment_correlation_matrix(df, segment_list)

        if args.strict:
            print("true")
            mask = (correlation_results_df >= cor_threshold).all(axis=1)
        else:
            print("false")
            mask = (correlation_results_df >= cor_threshold).any(axis=1)

        corr_df = correlation_results_df[mask]

        print(f"Write correlation matrix with a threshold of {cor_threshold}")
        corr_df.to_csv(args.output + ".tsv", sep="\t", index=True)
    else:
        df_transposed = df.transpose()
        correlation_matrix = create_correlation_matrix(df_transposed, args.p_value)

        print("Write correlation matrix.")
        correlation_matrix.to_csv(
            args.output + "_correlation_matrix.tsv", sep="\t", index=True
        )

        print("Write pairwise dataframe.")

        related_contigs = correlation_matrix[
            correlation_matrix >= cor_threshold
        ].stack()

        result_df = pd.DataFrame(related_contigs)
        result_df = result_df.reset_index()

        # Rename existing columns if necessary
        result_df.columns = ["Contig1", "Contig2", "Correlation"]
        result_df = result_df[result_df["Contig1"] != result_df["Contig2"]]

        result_df.sort_values(by="Contig1", inplace=True)
        result_df.to_csv(args.output + "_related_contigs.tsv", sep="\t", index=False)

    print("Finished.")


if __name__ == "__main__":
    main()
