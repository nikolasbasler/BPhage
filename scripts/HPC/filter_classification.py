import pandas as pd
import argparse

def get_last_non_empty_element(lst):
    for item in reversed(lst):
        if item.strip():
            return item

# Define command line arguments
parser = argparse.ArgumentParser(description='Filter genomad classification.')
parser.add_argument('genomad_virus_summary', type=str, help='genomad virus summary file path')
parser.add_argument('checkv_quality_summary', type=str, help='checkv quality summary file path')
parser.add_argument('ictv_vmr_xlxs', type=str, help='ICTV VMR Excel file path')
parser.add_argument('output_base_name', type=str, help='Output base name')
args = parser.parse_args()

# Load genomad output and append lowest classified taxon
bphage_ALL_1kb_cross_95_85_virus_summary = pd.read_csv(args.genomad_virus_summary, sep='\t')
bphage_ALL_1kb_cross_95_85_virus_summary['lowest_taxon'] = bphage_ALL_1kb_cross_95_85_virus_summary['taxonomy'].str.split(';').apply(get_last_non_empty_element)
# Load checkv output
quality_summary = pd.read_csv(args.checkv_quality_summary, sep='\t')

# Load ICTV taxonomy table
VMR_19_250422_MSL37 = pd.read_excel(args.ictv_vmr_xlxs)

# Make a list of phage-exclusive taxa
only_phages = VMR_19_250422_MSL37[(VMR_19_250422_MSL37['Host Source'] == "bacteria") | (VMR_19_250422_MSL37['Host Source'] == "archaea")]
excluding_phages = VMR_19_250422_MSL37[~((VMR_19_250422_MSL37['Host Source'] == "bacteria") | (VMR_19_250422_MSL37['Host Source'] == "archaea"))]
phage_exclusive_taxa = []
for level in ["Realm", "Kingdom", "Phylum", "Class", "Order", "Family"]:
    taxa = only_phages[level].dropna().unique()
    for tax in taxa:
        if tax not in excluding_phages[level].values:
            phage_exclusive_taxa.append(tax)

# Merge genomad and checkv outputs
merge_df = pd.merge(quality_summary, bphage_ALL_1kb_cross_95_85_virus_summary, left_on='contig_id', right_on='seq_name', how='inner')
merge_df = merge_df.drop(columns=['seq_name']) 

# Filter for Picobirnaviridae irrespective of CheckV's completeness estimate
picobirna_df = merge_df[merge_df['lowest_taxon'] == "Picobirnaviridae"]
picobirna_df = picobirna_df.fillna('NA')

# Filter merged output for at least 50% complete genomes
merge_filt_df = merge_df[merge_df['completeness'] >= 50]

# Filter for phage genomes
phages_df = merge_filt_df[merge_filt_df['lowest_taxon'].isin(phage_exclusive_taxa)]
phages_df = phages_df.fillna('NA')

# Filter for genomes classified only as virus or "Unclassified"
unclassified_viruses_df = merge_filt_df[(merge_filt_df['lowest_taxon'] == "Unclassified") | (merge_filt_df['lowest_taxon'] == "Viruses")]
unclassified_viruses_df = unclassified_viruses_df.fillna('NA')

# Write output csv files
picobirna_df.to_csv(args.output_base_name + "_picobirna.csv", index=False)
phages_df.to_csv(args.output_base_name + "_phages.csv", index=False)
unclassified_viruses_df.to_csv(args.output_base_name + "_unclassified_viruses.csv", index=False)


