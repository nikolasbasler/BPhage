library(ape)
library(vegan)
library(RLdbRDA)
source("scripts/R/custom_rldbrda.R")
library(tidyverse)

metadata <- readRDS("output/R/R_variables/metadata.RDS")

phage_tpm <- read.csv("output/R/relative_abundance/phage_tpm.csv") %>%
  tibble()

phold_predictions_with_extensions <- read.csv("output/R/gene_content/phold_predictions_with_extensions.csv") %>%
  tibble() %>%
  filter(str_starts(contig_id, "NODE"))

kegg_annotation <- read.delim("data/kegg_annotation.tsv", colClasses = "character") %>%
  tibble()

CDSs_with_metabolism_kegg <- kegg_annotation %>%
  # filter(Pathway_category == "Metabolism") %>%
  filter(Pathway_category %in% c("Metabolism", "Environmental Information Processing", NA)) %>%
  distinct(cds_id) %>% unlist(use.names = FALSE)


phold_predictions_with_extensions %>%
  left_join(kegg_annotation, ., by = "cds_id") %>% 
  select(cds_id, product, phrog, Pathway_category, Pathway_group, Pathway_name, Pathway_number, K_number) %>%
  filter(Pathway_category %in% c("Metabolism", "Environmental Information Processing", NA)) %>%
  View()

cropland_fraction <- read.csv("data/land_cover_results.csv") %>% 
  tibble() %>%
  mutate(cropland_fraction = cropland_fraction / 100) %>%
  rename(cropland_fraction_2k_radius = cropland_fraction) %>%
  arrange(cropland_fraction_2k_radius)

grene_presence_on_contigs <- phold_predictions_with_extensions %>% 
  filter(cds_id %in% CDSs_with_metabolism_kegg) %>%
  rename(contig = contig_id) %>%
  filter(function. == "moron, auxiliary metabolic gene and host takeover") %>%
  select(contig, product) %>% 
  mutate(present = 1) %>%
  distinct() %>% 
  mutate(product = str_replace_all(product, "-", "_"),
         product = str_replace_all(product, " ", "_")) %>%
  pivot_wider(names_from = product, values_from = present, values_fill = 0)

gene_tpm_long <- grene_presence_on_contigs %>%
  pivot_longer(-contig, names_to = "gene", values_to = "present") %>%
  left_join(., phage_tpm, by = "contig") %>%
  mutate(across(-c(contig, gene, present), ~ .x * present)) %>%
  select(-contig, -present) %>%
  pivot_longer(-gene, names_to = "Sample_ID", values_to = "tpm")

gene_presence <- list()
gene_presence$Sample_ID <- gene_tpm_long  %>%
  filter(tpm != 0) %>%
  group_by(gene, Sample_ID) %>%
  summarise(gene_presence = ifelse(sum(tpm) > 0 , 1, 0), .groups = "drop") %>%
  pivot_wider(id_cols = gene, names_from = Sample_ID, values_from = gene_presence, values_fill = 0)

for (thing in c("Bee_pool", "Country")) {
gene_presence[[thing]] <- gene_presence$Sample_ID %>%
  pivot_longer(-gene, names_to = "Sample_ID", values_to = "present") %>%
  left_join(metadata[c("Sample_ID", thing)], by = "Sample_ID") %>%
  group_by(gene, .data[[thing]]) %>%
  mutate(present = ifelse(sum(present) > 0, 1, 0)) %>%
  ungroup() %>%
  select(gene, all_of(thing), present) %>%
  distinct() %>%
  pivot_wider(id_cols = gene, names_from = all_of(thing), values_from = present)
}

gene_presence_gut_part <- list()
gene_presence_country_gut_part <- list()
for (g_part in c("mid", "ile", "rec")) {
  gene_presence_gut_part[[g_part]] <- gene_presence$Sample_ID %>%
    pivot_longer(-gene, names_to = "Sample_ID", values_to = "present") %>%
    filter(str_detect(Sample_ID, g_part)) %>%
    pivot_wider(id_cols = gene, names_from = Sample_ID, values_from = present)
  
  gene_presence_country_gut_part[[g_part]] <- gene_presence$Sample_ID %>%
    pivot_longer(-gene, names_to = "Sample_ID", values_to = "present") %>%
    filter(str_detect(Sample_ID, g_part)) %>%
    left_join(., metadata[c("Sample_ID", "Country")], by = "Sample_ID") %>%
    group_by(gene, Country) %>%
    mutate(present = ifelse(sum(present) > 0, 1, 0)) %>%
    ungroup() %>%
    select(gene, Country, present) %>%
    distinct() %>%
    pivot_wider(id_cols = gene, names_from = Country, values_from = present)
  }




# 
# gene_tpm <- list()
# gene_tpm$Sample_ID <- gene_tpm_long  %>%
#   filter(tpm != 0) %>%
#   group_by(gene, Sample_ID) %>%
#   summarise(gene_tpm = sum(tpm), .groups = "drop") %>%
#   pivot_wider(id_cols = gene, names_from = Sample_ID, values_from = gene_tpm, values_fill = 0)
# 
# for (thing in c("Bee_pool", "Country")) {
#   gene_tpm[[thing]] <- gene_tpm_long %>%
#     left_join(., metadata[c("Sample_ID", thing)], by = "Sample_ID") %>%
#     group_by(gene, .data[[thing]]) %>%
#     mutate(mean_tpm = mean(tpm)) %>%
#     ungroup() %>%
#     select(gene, all_of(thing), mean_tpm) %>%
#     distinct() %>%
#     pivot_wider(id_cols = gene, names_from = all_of(thing), values_from = mean_tpm)
# }

gene_presence$Bee_pool %>%
  pivot_longer(-gene) %>%
  group_by(name) %>%
  summarise(total_genes  = sum(value)) %>%
  arrange(total_genes)

# meta_filt <- gene_presence$Bee_pool %>%
#   pivot_longer(-gene, names_to = "Bee_pool") %>%
#   pivot_wider(id_cols = Bee_pool, names_from = gene, values_from = value) %>%
#   left_join(., metadata[c("Bee_pool", "Country", "Season")], by = "Bee_pool") %>%
#   select(Bee_pool, Country, Season, all_of(gene_presence$Bee_pool$gene)) %>%
#   distinct() %>%
#   as.data.frame()
# 
# gene_presence_dist <- gene_presence$Bee_pool %>%
#   select(-gene) %>%
#   t() %>%
#   vegdist(method = "jaccard")
# 
# ord <- pcoa(gene_presence_dist)
# 
# ord$vectors %>% 
#   as.data.frame() %>%
#   rownames_to_column("Bee_pool") %>%
#   as_tibble() %>%
#   select(Bee_pool, Axis.1, Axis.2) %>%
#   left_join(., metadata[c("Bee_pool", "Country")], by = "Bee_pool") %>%
#   ggplot(aes(x = Axis.1, y = Axis.2, color = Country)) +
#   geom_point()

gene_presence_dist <- gene_presence$Country %>% 
  select(-gene) %>%
  t() %>%
  vegdist(method = "jaccard")

meta_filt <- gene_presence$Country %>%
  pivot_longer(-gene) %>%
  pivot_wider(id_cols = name, names_from = gene, values_from = value) %>%
  rename(Country = name) %>%
  distinct() %>%
  as.data.frame()

ord <- pcoa(gene_presence_dist)

ord$vectors %>% 
  as.data.frame() %>%
  rownames_to_column("Country") %>%
  as_tibble() %>%
  select(Country, Axis.1, Axis.2) %>%
  left_join(., metadata[c("Bee_pool", "Country")], by = "Country") %>%
  ggplot(aes(x = Axis.1, y = Axis.2, color = Country)) +
  geom_point()

  

RDAs <- custom_rldbrda(gene_presence_dist, meta_filt)

plot_data <- RLdbRDA::prepare_plot_data(RDAs)

RDA_plot <- RLdbRDA::plot_dbrda(plot_data) + 
  scale_fill_manual(values=c("#ef8f01", "#8B4513"),
                    labels=c(bquote(R^2), bquote('Cumulative' ~ R^2))) +
  # scale_x_continuous(limits = c(-10, 10)) +
  # ggtitle(paste0(tax, "_", plotting_lable[[set]])) +
  theme(legend.position = "bottom")
