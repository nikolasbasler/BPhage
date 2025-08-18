# Please run renv::restore() in order to install all required packages in the 
# versions used here.

library(phyloseq)
library(decontam)
library(tidyverse)
library(vegan)

source("scripts/R/helpers/tpm_and_reads_per_kb.R")
source("scripts/R/helpers/decontam.R")

abundance.table <- read.csv("output/mapping_stats_phages/stats.phages.mapped_reads.csv", row.names=1) %>%
  rownames_to_column("contig") %>%
  filter(contig != "NODE_A1975_length_2506_cov_68.193907_PT_19410_aut_rec_d") %>% # Filter out the one Picobirna contig that isn't the RdRp segment
  column_to_rownames("contig")
horizontal.cov.table <- read.csv("output/mapping_stats_phages/stats.phages.horizontal_coverage.csv", row.names=1) %>%
  rownames_to_column("contig") %>%
  filter(contig != "NODE_A1975_length_2506_cov_68.193907_PT_19410_aut_rec_d") %>% # Filter out the one Picobirna contig that isn't the RdRp segment
  column_to_rownames("contig")
mean_depth_table <- read.csv("output/mapping_stats_phages/stats.phages.mean_depth.csv", row.names=1) %>%
  rownames_to_column("contig") %>%
  filter(contig != "NODE_A1975_length_2506_cov_68.193907_PT_19410_aut_rec_d") %>% # Filter out the one Picobirna contig that isn't the RdRp segment
  column_to_rownames("contig")

gnmd_classification <- read.csv("output/bphage_ALL_1kb_phages.csv") %>%
  rbind(read.csv("output/bphage_ALL_1kb_unclassified_viruses.csv")) %>%
  rbind(read.csv("output/bphage_ALL_1kb_picobirna.csv")) %>% 
  separate_wider_delim(taxonomy, delim = ";",
                       names = c("Classification", "Realm", "Kingdom", "Phylum",
                                 "Class", "Order", "Family"),
                       too_few = "align_start") %>%
  mutate(Genus = "Unclassified",
         Species= "Unclassified", .before = lowest_taxon) %>%
  mutate(across(c(Realm, Kingdom, Phylum, Class, Order, Family), ~ ifelse(is.na(.x) | .x == "", "Unclassified", .x))) %>%
  rename(contig=contig_id) %>%
  filter(contig != "NODE_A1975_length_2506_cov_68.193907_PT_19410_aut_rec_d") # %>% # Filter out the one Picobirna contig that isn't the RdRp segment

metadata <- readRDS("data/metadata.RDS") %>% as.data.frame()
row.names(metadata) <- metadata$Sample_ID # decontam needs this...

contig_length_df <- read.csv("output/mapping_stats_phages/phages.contig_lentghs.csv") %>%
  mutate(length = length/1000) %>%
  rename(length_kb = length) %>%
  filter(contig != "NODE_A1975_length_2506_cov_68.193907_PT_19410_aut_rec_d") # Filter out the one Picobirna contig that isn't the RdRp segment

################################################################################
################################################################################
# Filtering by horizontal coverage and mean depth

hzc.filter=70
mean.depth.filter = 1

if (!identical(rownames(abundance.table), rownames(horizontal.cov.table)) |
    !identical(colnames(abundance.table), colnames(horizontal.cov.table)) |
    !identical(rownames(abundance.table), rownames(mean_depth_table)) |
    !identical(colnames(abundance.table), colnames(mean_depth_table))) 
{
  stop("Row or column names in hzc or mean depth table don't match up!")
}

abundance_table_filt <- abundance.table
abundance_table_filt[horizontal.cov.table < hzc.filter] <- 0
abundance_table_filt[mean_depth_table < mean.depth.filter] <- 0

################################################################################
# Filter out contaminants
decontam <- decontam_identification(abtable=abundance_table_filt, 
                                    classification=gnmd_classification)

# NODE_A1_length_19355_cov_19.948802_Blank_pool_01 found in all blank pools
# (except WTA blanks) and in 73 samples
# NODE_A1_length_30014_cov_53.072886_Blank_pool_12 found in 3 blank pools
# and 1 sample.

abundance_table_filt <- abundance_table_filt %>%
  rownames_to_column("contig") %>%
  filter(!contig %in% decontam$contaminants) %>%
  column_to_rownames("contig") 

################################################################################
# Filter out contigs that, after the hzc filter, are only present in blanks,
# then filter out contigs and samples that only have zeros.
contigs_before_hzc <- rownames(abundance_table_filt)
samples_before_hzc <- colnames(abundance_table_filt)

contigs_only_in_blanks <- abundance_table_filt %>%
  select(!contains("Blank")) %>%
  filter(rowSums(.)==0) %>%
  rownames()
abundance_table_filt <- abundance_table_filt %>%
  rownames_to_column("contig") %>%
  filter(!contig %in% contigs_only_in_blanks) %>%
  column_to_rownames("contig") %>%
  select(!names(.)[colSums(.) == 0]) %>%
  filter(rowSums(.)>0) 

contigs_after_hzc = rownames(abundance_table_filt)
removed_contigs = setdiff(contigs_before_hzc, contigs_after_hzc)

samples_after_hzc = colnames(abundance_table_filt)
removed_samples = setdiff(samples_before_hzc, samples_after_hzc)

remaining_pools <- metadata %>%
  mutate(pool = paste(Country, Hive_ID, Season, sep="_")) %>%
  select(Sample_ID, pool) %>%
  filter(!Sample_ID %in% removed_samples) %>%
  select(pool) %>%
  distinct() %>%
  nrow()

cat(length(contigs_before_hzc)-length(contigs_after_hzc), "of", length(contigs_before_hzc),
          "contigs and", length(samples_before_hzc)-length(samples_after_hzc), "of", length(samples_before_hzc),
          "samples removed due to horizontal coverage filter of", hzc.filter,", mean depth filter of", mean.depth.filter, "or because only present in blanks.\n",
          "\nRemoved contigs:\n", removed_contigs, "\n",
          "\nRemoved samples:\n", removed_samples, "\n",
          "\nRemaining bee pools:", remaining_pools)

cat("Lowest number of mapped reads (>0) to any contig:", min(abundance_table_filt[abundance_table_filt>0]),"\n")

################################################################################
################################################################################

# Write output files
system("mkdir -p output/R/")

# Decontam graph
ggsave("output/R/decontam.library.size.by.control.or.sample.pdf", decontam$plot, width=30, height=10)

# Raw abundance tables
write_csv(rownames_to_column(abundance_table_filt, "contig"), file = "output/R/phage.filt.abundance.contig.csv")

# Classification table
gnmd_classification %>%
  filter(contig %in% rownames(abundance_table_filt)) %>%
  write_csv(file = "output/R/phage.filt.gnmd.classification.csv")



