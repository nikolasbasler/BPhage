library(tidyverse)
library(patchwork)
library(ggVennDiagram)

metadata <- readRDS("data/metadata.RDS") %>% tibble()
classification <- readRDS("data/classification.RDS")

phage_tpm <- read.csv("output/R/relative_abundance/phage_tpm.csv") %>%
  tibble()
present_in_all_countries <- read_lines("data/core_contigs.txt")

map_stats <- read.delim("output/bphage_viper_output/read_stats.tsv") %>%
  tibble()

filt_abundance_table <- read.csv("output/R/phage.filt.abundance.contig.csv") %>%
  tibble()
abundance_table <- read.csv("output/mapping_stats_phages/stats.phages.mapped_reads.csv") %>% 
  tibble()

bphage_and_inpha_95.85_clusters <- read.delim("output/inphared_clustering/bphage_and_inpha_70-85_clusters.tsv", header=FALSE) %>%
  tibble()

##### 
## Calculating this factor of 25, regarding the value bee pollination services vs. the value of honey and beeswax:

# according to the federal reserve bank of st. louis (https://fred.stlouisfed.org/graph/?g=RJJP)
exchange_rate_2005 <- read.csv("data/Dollar_to_Euro_exchange_rates_2005_from_fred.stlouisfed.org.csv") %>%
  filter(!is.na(DEXUSEU)) %>%
  reframe(mean_rate = mean(DEXUSEU)) %>%
  unlist(use.names = FALSE)

pollination_value_euro <- 153e9/2 # Gallai et al 2009 estimate 153 billion euro by insect pollination in 2005. Reilly et al 2024 estimate honey bees contribute about half.
pollination_value_dollar <- pollination_value_euro * exchange_rate_2005

FAOSTAT_bee_product_value_en_5.30.2025 <- read.csv("data/FAOSTAT_bee_product_value_en_5-30-2025.csv") %>%
  tibble()

bee_product_value <- FAOSTAT_bee_product_value_en_5.30.2025 %>%
  select(Item, Value) %>%
  reframe(bee_product_value = sum(Value)*1000) %>%
  unlist(use.names = FALSE)

pollination_value_dollar / bee_product_value

##### 
## Venn diagram of core phages shared in different gut parts
presence <- phage_tpm %>% 
  pivot_longer(-contig, names_to = "Sample_ID", values_to = "tpm") %>%
  filter(contig %in% present_in_all_countries) %>%
  mutate(present = tpm > 0) %>%
  left_join(., metadata[c("Sample_ID", "Gut_part")], by = "Sample_ID")

for_venn <- list()
for (gpart in unique(presence$Gut_part)) {
  for_venn[[gpart]] <- presence %>% 
    filter(Gut_part == gpart,
           present) %>%
    distinct(contig) %>%
    unlist(use.names = FALSE)
}

core_shared_between_guts <- ggVennDiagram(for_venn)

##### 
## total observed contigs per bee pool
contigs_per_bee_pool <- phage_tpm %>%
  pivot_longer(-contig, names_to = "Sample_ID", values_to = "tpm") %>%
  left_join(., metadata[c("Sample_ID", "Bee_pool")], by = "Sample_ID") %>%
  group_by(Bee_pool) %>%
  summarise(contig_count = sum(tpm > 0))

contigs_per_bee_pool_hist <- ggplot(contigs_per_bee_pool, aes(x = contig_count)) +
  geom_histogram(bins = 50) +
  geom_vline(xintercept = mean(contigs_per_bee_pool$contig_count)) +
  labs(x = "Phage contigs per bee pool")

##### 
## MAPPING STATS SUMMARY:
pre_mapping_stats <- map_stats %>%
  rename("Sample_ID" = Sample) %>%
  filter(!str_detect(Sample_ID, "Blank")) %>%
  pivot_longer(-Sample_ID, names_to = "metric", values_to = "raw_reads") %>%
  mutate(metric = factor(metric, levels = c("Raw_R1_plus_R2", "Deduplicated_R1_plus_R2", 
                                            "Trimmed_R1_plus_R2", "Trimmed_unpaired", "Trimmed_total",
                                            "Hostout_R1_plus_R2", "Hostout_unpaired", "Hostout_total"))) %>%
  group_by(metric) %>%
  summarise(total = sum(raw_reads),
            average = mean(raw_reads),
            total_in_millions = total / 1000000,
            average_in_millions = average / 1000000) 

phage.filt.abundance.contig_long <- filt_abundance_table %>%
  pivot_longer(-contig, names_to = "Sample_ID", values_to = "filtered_mapped")
abundance.table_long <- abundance_table %>% 
  filter(contig != "NODE_A1975_length_2506_cov_68.193907_PT_19410_aut_rec_d") %>% # Filter out the one Picobirna contig that isn't the RdRp segment
  pivot_longer(-contig, names_to = "Sample_ID", values_to = "unfilterd_mapped") %>%
  filter(!str_detect(Sample_ID, "Blank"))

post_mapping_stats <- full_join(phage.filt.abundance.contig_long, abundance.table_long, by = c("contig", "Sample_ID")) %>%
  group_by(Sample_ID) %>%
  summarise(unfilterd_mapped = sum(unfilterd_mapped),
            filtered_mapped = sum(filtered_mapped, na.rm = TRUE)) %>%
  reframe(total_unfiltered = sum(unfilterd_mapped),
          avaerage_unfiltered = mean(unfilterd_mapped),
          total_filtered = sum(filtered_mapped, na.rm = TRUE),
          average_filtered = mean(filtered_mapped, na.rm = TRUE)) %>%
  pivot_longer(everything(), names_to = "metric", values_to = "mapped_reads") %>%
  mutate(mapped_reads_in_millions = mapped_reads/1000000)

pre_mapping_stats
post_mapping_stats

##### 
## INPHARED CLUSTERING

relevant_clusters_tibble <- bphage_and_inpha_95.85_clusters %>%
  filter(str_detect(V2, "NODE"),
         str_detect(V2, ",")) %>%
  rename("representative" = V1,
         "members" = V2)

relevant_clusters_list <- str_split(relevant_clusters_tibble$members, ","); names(relevant_clusters_list) <- relevant_clusters_tibble$representative

clustered_nodes <- relevant_clusters_list %>%
  lapply(function(x) grep("NODE", x, value = TRUE)) %>%
  unlist(use.names = FALSE)

contigs_with_inphas_list <- list()
for (node in clustered_nodes) {
    temp <- relevant_clusters_tibble %>%
    filter(str_detect(members, node)) %>%
    mutate(members = gsub("NODE_[^,]*", "", members),
           members = gsub(",{2,}", ",", members),
           members = gsub("^,|,$", "", members)
           ) %>%
    select(members) %>%
    unlist() %>%
    str_split(",") %>%
    set_names(node)
    
    if (all(temp[[node]] != "")) {
      contigs_with_inphas_list <- c(contigs_with_inphas_list, temp)
    }
}

inpha_clustering <- tibble(contig = names(contigs_with_inphas_list),
                                     INPHARED_clustered = lapply(contigs_with_inphas_list, function(x) paste(x, collapse = ",")) %>%
                                       unlist())

new_classification_df <- classification %>%
  select(!starts_with("INPHARED")) %>%
  left_join(., inpha_clustering, by = "contig") %>% 
  relocate(INPHARED_clustered, .before = Prevalence_other_datasets) %>%
  rename(INPHARED_and_Bueren_clustered = INPHARED_clustered)

#####
#### Save files

# Venn diagram
ggsave("output/R/core_shared_between_guts.pdf", core_shared_between_guts,
       width = 5, height = 5)

# Histogram contigs per bee pool
write_csv(contigs_per_bee_pool, "output/R/alpha/total_contigs_per_bee_pool.csv")
ggsave("output/R/alpha/total_contigs_per_bee_pool.pdf", contigs_per_bee_pool_hist,
       width = 5, height = 5)

# For convenience, to avoid backtracking
# write_csv(new_classification_df, "output/R/classification.csv")
