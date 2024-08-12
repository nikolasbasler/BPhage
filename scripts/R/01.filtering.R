library(phyloseq)
library(decontam)
library(tidyverse)
library(vegan)

source("scripts/R/helpers/tpm_and_reads_per_kb.R")
source("scripts/R/helpers/decontam.R")

# THIS PATH HAS TO BE ADJUSTED BY THE USER! PANDAS AND SCIPY HAVE TO BE INSTALLED
# This is only for the co-occurrence analysis, though. Might be skipped anyway
python_path="/Users/nikolasbasler/miniforge3/envs/co_occurrence/bin/python3"

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

# Abandoned apporach:
# abundance.table <- read.csv("output/mapping_stats_bphage_and_others/stats.bphage_and_others.mapped_reads.csv", row.names=1) %>%
#   rownames_to_column("contig") %>%
#   filter(contig != "NODE_A1975_length_2506_cov_68.193907_PT_19410_aut_rec_d") %>% # Filter out the one Picobirna contig that isn't the RdRp segment
#   column_to_rownames("contig")
# horizontal.cov.table <- read.csv("output/mapping_stats_bphage_and_others/stats.bphage_and_others.horizontal_coverage.csv", row.names=1) %>%
#   rownames_to_column("contig") %>%
#   filter(contig != "NODE_A1975_length_2506_cov_68.193907_PT_19410_aut_rec_d") %>% # Filter out the one Picobirna contig that isn't the RdRp segment
#   column_to_rownames("contig")
# mean_depth_table <- read.csv("output/mapping_stats_bphage_and_others/stats.bphage_and_others.mean_depth.csv", row.names=1) %>%
#   rownames_to_column("contig") %>%
#   filter(contig != "NODE_A1975_length_2506_cov_68.193907_PT_19410_aut_rec_d") %>% # Filter out the one Picobirna contig that isn't the RdRp segment
#   column_to_rownames("contig")

gnmd_classification <- read.csv("output/bphage_ALL_1kb_phages.csv") %>%
  rbind(read.csv("output/bphage_ALL_1kb_unclassified_viruses.csv")) %>%
  rbind(read.csv("output/bphage_ALL_1kb_picobirna.csv")) %>% 
  # rbind(read.csv("output/other_studies_phages.csv")) %>%
  # rbind(read.csv("output/other_studies_unclassified_viruses.csv")) %>%
  separate_wider_delim(taxonomy, delim = ";",
                       names = c("Classification", "Realm", "Kingdom", "Phylum",
                                 "Class", "Order", "Family"),
                       too_few = "align_start") %>%
  mutate(Genus = "Unclassified",
         Species= "Unclassified", .before = lowest_taxon) %>%
  mutate(across(c(Realm, Kingdom, Phylum, Class, Order, Family), ~ ifelse(is.na(.x) | .x == "", "Unclassified", .x))) %>%
  rename(contig=contig_id) %>%
  filter(contig != "NODE_A1975_length_2506_cov_68.193907_PT_19410_aut_rec_d") # %>% # Filter out the one Picobirna contig that isn't the RdRp segment
  # filter(contig %in% rownames(abundance.table))

host_groups_df<- read_csv("data/host_groups.csv", show_col_types = FALSE)

metadata=metadata <- read.csv("data/metadata.csv")
row.names(metadata) <- metadata$Sample_ID

contig_length_df <- read.csv("output/mapping_stats_phages/phages.contig_lentghs.csv") %>%
  mutate(length = length/1000) %>%
  rename(length_kb = length) %>%
  filter(contig != "NODE_A1975_length_2506_cov_68.193907_PT_19410_aut_rec_d") # Filter out the one Picobirna contig that isn't the RdRp segment

################################################################################
################################################################################
## Add host group information to classification table

gnmd_classification <- host_groups_df %>%
  select(Rank, Taxon, Host_group) %>%
  left_join(gnmd_classification, ., by=join_by(Classification==Taxon)) %>%
  select(-Rank) %>%
  mutate(Host_group = ifelse(is.na(Host_group), "unassigned", Host_group))


# contig_length_df %>%
#   ggplot(aes(x=length_kb)) +
#   geom_histogram(binwidth=0.5) +
#   geom_vline(xintercept = 3, color="red")

# Min contig length filter
# min_kb_filter = 3
# abundance.table <- contig_length_df %>%
#   filter(length_kb >= min_kb_filter) %>%
#   inner_join(., rownames_to_column(abundance.table, "contig"), by="contig") %>% 
#   select(-length_kb) %>%
#   column_to_rownames("contig")

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
## Co-occurrence analysis

system("mkdir -p output/R/co_occurrence/")
abundance_table_filt %>%
  rownames_to_column("contig") %>% 
  select(contig, contains("_rec_")) %>% 
  arrange(contig) %>%
  write_tsv("output/R/co_occurrence/abundances.for.co-occ.tsv")

contig_length_df %>%
  filter(contig %in% rownames(abundance_table_filt)) %>%
  arrange(contig) %>%
  write_tsv("output/R/co_occurrence/contig.lengths.tsv", col_names = FALSE)

system(paste0(python_path,
              " scripts/R/helpers/co-occurrence.py",
              " --input output/R/co_occurrence/abundances.for.co-occ.tsv",
              " --lengths output/R/co_occurrence/contig.lengths.tsv",
              " --output output/R/co_occurrence/co_occurrences",
              " --correlation 0.5"))
co_occurrences_related_contigs <- read.delim("output/R/co_occurrence/co_occurrences_related_contigs.tsv")

co_occ_join <- co_occurrences_related_contigs %>%
  inner_join(., gnmd_classification, by=join_by(Contig1==contig)) %>%
  inner_join(., gnmd_classification, by=join_by(Contig2==contig)) %>%
  select(-contains(c("Avg_log_e.value", "Classification", 
                     "Classified_by"))) %>%
  select(Contig1, Contig2, Correlation, completeness.x, topology.x, virus_score.x, n_hallmarks.x, lowest_taxon.x,
         completeness.y, topology.y, virus_score.y, n_hallmarks.y, lowest_taxon.y) %>%
  rename(completeness_1 = completeness.x,
         topology_1 = topology.x,
         virus_score_1 = virus_score.x,
         n_hallmarks_1 = n_hallmarks.x,
         lowest_taxon_1 = lowest_taxon.x,
         completeness_2 = completeness.y,
         topology_2 = topology.y,
         virus_score_2 = virus_score.y,
         n_hallmarks_2 = n_hallmarks.y,
         lowest_taxon_2 = lowest_taxon.y)

co_occ_stats <- co_occ_join %>%
  select(lowest_taxon_1, lowest_taxon_2) %>%
  group_by_all() %>%
  mutate(number_of_correlations=n()) %>% 
  unique() %>%
  arrange(desc(number_of_correlations))

################################################################################
################################################################################
# TPM and prevalence filter (Not in effect right now)

contig_tpm_temp <- rownames_to_column(abundance_table_filt, "contig") %>%
  calc_tpm(., "contig", contig_length_df)
  
gnmd_classification_refined <- gnmd_classification # For compatiblity further down

# tpm_thresh <- 0.01
# prev_thresh <- 1
# 
# passed_tpm_and_prevalence_filters <- contig_tpm_temp %>%
#   select(!contains("Blank")) %>%
#   column_to_rownames("contig") %>%
#   filter(if_any(everything(), ~ . >= tpm_thresh)) %>%
#   filter(rowSums(. > 0) >= prev_thresh) %>%
#   rownames_to_column("contig") %>%
#   inner_join(., gnmd_classification, by="contig") %>%
#   select(contig) %>%
#   unlist(use.names = FALSE)
# 
# gnmd_classification_refined <- gnmd_classification %>%
#   filter(contig %in% passed_tpm_and_prevalence_filters) %>%
#   column_to_rownames("contig")
# 
# abundance_table_filt <- rownames_to_column(abundance_table_filt, "contig") %>%
#   filter(contig %in% passed_tpm_and_prevalence_filters) %>%
#   column_to_rownames("contig")

##
max_tpm_and_prev_df <- contig_tpm_temp %>%
  select(!contains("Blank")) %>%
  rowwise() %>%
  mutate(max_tpm = max(c_across(-contig)),
         mean_tpm_where_present = rowMeans(across(-c(contig, max_tpm), ~ifelse(. > 0, ., NA)), na.rm = TRUE),
         prevalence = sum(c_across(-c(contig, max_tpm, mean_tpm_where_present)) > 0),
         # prevalence_of_tpm_5perc = sum(c_across(-c(contig, max_tpm, mean_tpm_where_present, prevalence)) > 0.05),
         prevalence_of_tpm_1perc = sum(c_across(-c(contig, max_tpm, mean_tpm_where_present, prevalence)) > 0.01)) %>%
  select(contig, max_tpm, mean_tpm_where_present, prevalence, prevalence_of_tpm_1perc) %>%
  inner_join(gnmd_classification_refined,. , by="contig")

max_tpm_and_prev_hist <- max_tpm_and_prev_df %>%
  select(contig, max_tpm, mean_tpm_where_present, prevalence, prevalence_of_tpm_1perc) %>%
  pivot_longer(-contig) %>%
  ggplot(aes(x=value)) +
  geom_histogram(bins=50) +
  facet_wrap(~name, scales="free")

################################################################################
################################################################################

# Write output files

# Decontam graph
ggsave("output/R/decontam.library.size.by.control.or.sample.pdf", decontam$plot, width=30, height=10)

# Co-occurrence
write.csv(co_occ_join, "output/R/co_occurrence/co_occ_joined.csv", row.names=FALSE)
write.csv(co_occ_stats, "output/R/co_occurrence/co_occ_stats.csv", row.names=FALSE)

# Max TPM and prevalance
write.csv(max_tpm_and_prev_df, "output/R/max_tpm_and_prev.csv", row.names=FALSE)
ggsave("output/R/max_tpm_and_prev.pdf", max_tpm_and_prev_hist, width=8, height=5)

# Raw abundance tables
write_csv(rownames_to_column(abundance_table_filt, "contig"), file = "output/R/phage.filt.abundance.contig.csv")

# Classification table
gnmd_classification_refined %>%
  filter(contig %in% rownames(abundance_table_filt)) %>%
  write_csv(file = "output/R/phage.filt.gnmd.classification.csv")



