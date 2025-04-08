library(tidyverse)
library(ggpubr)
library(patchwork)

metadata <- readRDS("output/R/R_variables/metadata.RDS")
classification <- readRDS("output/R/R_variables/classification.RDS")
present_in_all_countries <- read_lines("data/core_contigs.txt")

# all_hosts <- read.csv("output/R/host_pies/all_hosts.all.csv") %>%
#   rename(contig = Virus, Host_genus = Genus)

phold_predictions_with_extensions <- read.csv("output/R/gene_content/phold_predictions_with_extensions.csv") %>%
  tibble() %>%
  filter(str_starts(contig_id, "NODE"))

cropland_fraction <- read.csv("data/land_cover_results.csv") %>% 
  tibble() %>%
  mutate(cropland_fraction = cropland_fraction / 100) %>%
  rename(cropland_fraction_2k_radius = cropland_fraction) %>%
  arrange(cropland_fraction_2k_radius)

absolute_counts <- read.csv("output/R/absolute_counts.csv") %>% 
  tibble()


moron_gene_copies <- phold_predictions_with_extensions %>%
  rename("contig" = "contig_id") %>%
  filter(str_detect(contig, "_rec")) %>% 
  filter(function. == "moron, auxiliary metabolic gene and host takeover") %>%
  select(contig, cds_id, function., product) %>%
  left_join(., absolute_counts, by = "contig") %>%
  select(-contig) %>%
  filter(if_all(-cds_id, ~ !is.na(.))) %>%
  pivot_longer(-c(cds_id, product, function.), names_to = "Sample_ID", values_to = "VLPs_per_ul") %>%
  # left_join(., metadata[c("Sample_ID", "Country")], by = "Sample_ID") %>%
  # group_by(product, Country) %>%  
  left_join(., metadata[c("Sample_ID", "Hive_ID")], by = "Sample_ID") %>%
  group_by(product, Hive_ID) %>%
  summarise(gene_copies_per_ul = sum(VLPs_per_ul), .groups = "drop") %>%
  left_join(., metadata[c("Hive_ID", "Country")] %>% distinct(), by = "Hive_ID")

gene_cor <- list()
for (gene in unique(moron_gene_copies$product)) {
  gene_cor[[gene]] <- moron_gene_copies %>%
    filter(product == gene) %>%
    left_join(., cropland_fraction, by = "Country") %>%
    ggplot(aes(x = cropland_fraction_2k_radius, y = gene_copies_per_ul)) +
    geom_point() +
    geom_smooth(method = "glm", formula = y ~ x) +
    stat_cor(method = "spearman") +
    ggtitle(gene)
}

