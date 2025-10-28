library(tidyverse)
library(ggpubr)
library(ggtext)
library(patchwork)
library(forcats)

source("scripts/R/helpers/mixed_helpers.R")

metadata <- readRDS("data/metadata.RDS")
classification <- readRDS("data/classification.RDS")
present_in_all_countries <- read_lines("data/core_contigs.txt")
phage_tpm <- read.csv("output/R/relative_abundance/phage_tpm.csv") %>%
  tibble()

phold_annotations_extended <- read.delim("output/core_contig_refinement/extended_contigs_phold/phold_per_cds_predictions_long_names.tsv") %>%
  tibble() %>%
  mutate(contig_id = str_extract(contig_id, "^([^_]+_){10}[^_]+")) %>%
  mutate(cds_id = paste0(contig_id, "_", str_extract(cds_id, "[^_]+_[^_]+$")))

phold_annotations_unextended <- read.delim("output/annotation/phold_compare_bphage_and_others/phold_per_cds_predictions_long_names.tsv") %>%
  tibble()

## TODO: add these files to the midsave
vfdb_cds_predictions_extended <- read.delim("output/core_contig_refinement/extended_contigs_phold/sub_db_tophits/vfdb_cds_predictions_long_names.tsv") %>%
  tibble()
vfdb_cds_predictions_unextended <- read.delim("output/annotation/phold_compare_bphage_and_others/sub_db_tophits/vfdb_cds_predictions_long_names.tsv") %>%
  tibble() %>%
  filter(!contig_id %in% vfdb_cds_predictions_extended$contig_id)
vfdb_cds_predictions_with_extensions <- rbind(vfdb_cds_predictions_extended, vfdb_cds_predictions_unextended)

## TODO: if this is to be loaded here, change order of scripts!
# AMG curation would have to be run first then

## CONTINUE HERE: UPDATE THIS TABLE SO IT ALSO CONTAINS contig I!
remove_for_stringency_PAPS <- read.delim("output/R/AMG_curation/remove_for_stringency.PAPS reductase.tsv")


#

phold_predictions_with_extensions <- phold_annotations_unextended %>%
  filter(str_starts(contig_id, "NODE"),
         !contig_id %in% phold_annotations_extended$contig_id) %>%
  rbind(., phold_annotations_extended) %>%
  # filter(str_starts(contig_id, "NODE")) %>%
  filter(contig_id %in% classification$contig) %>%
  mutate( # Some of this is pre-peer review legacy:
    product = str_replace_all(product, "levanase", "Levanase"),
    product = str_replace_all(product, "glutamine amidotransferase", "GATase"),
    product = str_replace_all(product, "glucosyltransferase", "Glucosyltransferase"),
    product = str_replace_all(product, "porphyrin biosynthesis", "Porphyrin synthesis protein"),
    product = str_replace_all(product, "aerobic cobaltochelatase CobT subunit", "CobT"),
    product = str_replace_all(product, "chitinase", "Chitinase"),
    product = str_replace_all(product, "dTDP-4-dehydrorhamnose 3", "RmlC"),
    product = str_replace_all(product, "nicotinamide-nucleotide adenylyltransferase", "NMNAT"),
    product = str_replace_all(product, "PnuC-like nicotinamide mononucleotide transport", "PnuC"),
    product = str_replace_all(product, "phosphoadenosine phosphosulfate reductase", "PAPS reductase"),
    product = str_replace_all(product, "NrdD-like anaerobic ribonucleotide reductase large subunit", "NrdD"),
    product = str_replace_all(product, "ribosomal protein S6 glutaminyl transferase", "RimK")
    ) %>%
  left_join(., vfdb_cds_predictions_with_extensions[c("cds_id", "description")], by = "cds_id") %>%
  mutate(product = ifelse(is.na(description), product, paste0("VFDB ", description))) %>%
  select(-description)

kegg_mapping <- read.delim("data/kegg_mapping.tsv", colClasses = "character") %>%
  tibble()

kegg_and_phold <- kegg_mapping %>%
  left_join(., phold_predictions_with_extensions[c("contig_id", "cds_id", "phrog", "function.", "product")], by = "cds_id")

phage_tpm <- read.csv("output/R/relative_abundance/phage_tpm.csv") %>%
  tibble()

vibrant_amg_KOs <- read_lines("data/AMG_vibrant_curated_KOs")

kegg_and_phold %>%
  summarise(cds_with_k = n_distinct(cds_id))

kegg_and_phold %>%
  filter(K_number %in% vibrant_amg_KOs) %>%
  summarise(cds_with_k = n_distinct(cds_id))


#####
# PHROG BARS
# moron_bar_colors <- rev(c("lightgrey", "#E4B3E4", "#555555"))
moron_bar_colors <- rev(c("lightgrey", "#F9DDF9", "#555555"))

phrog_tibble <- phold_predictions_with_extensions %>%
  rename(contig = contig_id) %>%
  left_join(., classification[c("contig", "Core")], by = "contig") %>%
  group_by(function., Core) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(function.) %>%
  mutate(count_all = sum(count)) %>%
  ungroup() %>%
  filter(Core == "yes") %>%
  select(-Core) %>%
  rename(count_core = count) %>%
  mutate(collapsed_cat = case_when(function. == "unknown function" ~"unknown function",
                                   function. == "moron, auxiliary metabolic gene and host takeover" ~ '"moron", AMG\nand host takeover',
                                   .default = "other PHROG\ncategories"
                                   )) %>%
  group_by(collapsed_cat) %>%
  summarise(core = sum(count_core),
            all = sum(count_all)) %>%
  mutate(collapsed_cat = factor(collapsed_cat, levels = c("other PHROG\ncategories", '"moron", AMG\nand host takeover', "unknown function")))

phrog_bar_vertical <- list()
phrog_bar_vertical$both <- phrog_tibble %>%
  pivot_longer(-collapsed_cat) %>%
  mutate(name = factor(name, levels = c("core", "all"))) %>%
  ggplot(aes(x = name, y = value, fill = collapsed_cat)) +
  geom_col() +
  theme_void() +
  theme(
    axis.text.x  = element_text(vjust = 5),
    legend.margin=margin(0,0,0,10),
    legend.position = "left",
    legend.key.spacing.y = unit(3, "pt") 
    ) +
  scale_fill_manual(values = moron_bar_colors) +
  labs(fill = "PHROG category") +
  guides(fill = guide_legend(label.position = "left"))

for (set in c("all", "core")) {
  phrog_bar_vertical[[set]] <- phrog_tibble %>%
    pivot_longer(-collapsed_cat) %>%
    filter(name == set) %>%
    ggplot(aes(
      x    = name,
      y    = value,
      fill = collapsed_cat
    )) +
    geom_col(position = "fill") +
    scale_fill_manual(values = moron_bar_colors) +
    theme_void() +
    theme(
      legend.position = "none",
      plot.margin = margin(5, 5, 5, 5, "pt")
    )
}

phrog_bar_horizontal <- list()
for(set in names(phrog_bar_vertical)) {
  phrog_bar_horizontal[[set]] <- phrog_bar_vertical[[set]] + 
    coord_flip() +
    theme(legend.position = "top",
          legend.text = element_text(hjust = 0),
          plot.margin = margin(t = 20, unit = "pt"),
          axis.text.y  = element_text()
          ) +
    guides(fill = guide_legend(reverse=T))
}

#####
# KEGG BAR

CDSs_with_metabolism_kegg <- kegg_and_phold %>%
  filter(K_number %in% vibrant_amg_KOs) %>%
  select(cds_id, contig_id, product) %>%
  distinct()

genes_with_kegg <- phold_predictions_with_extensions %>% 
  filter(cds_id %in% CDSs_with_metabolism_kegg$cds_id) %>%
  distinct(product) %>%
  unlist(use.names = FALSE)

kegg_tibble <- phold_predictions_with_extensions %>%
  filter(product %in% genes_with_kegg) %>%
  reframe(placeholder = "placeholder",
          `Assigned metabol. gene` = length(CDSs_with_metabolism_kegg$cds_id),
          `other assigned` = n_distinct(kegg_mapping$cds_id) - `Assigned metabol. gene`,
          `Unass. metabol. gene` = n_distinct(cds_id) - `Assigned metabol. gene`,
          `other unassigned` = phold_predictions_with_extensions %>% 
            filter(function. == "moron, auxiliary metabolic gene and host takeover") %>% 
            n_distinct("cds_id") - sum(`Assigned metabol. gene`, `other assigned`, `Unass. metabol. gene`)
          ) %>%
  pivot_longer(-placeholder, names_to = "mapping", values_to = "gene_count") %>%
  select(-placeholder) %>%
  mutate(mapping_label = paste0(mapping, " (", gene_count, ")")) %>%
  mutate(mapping = factor(mapping, levels = c("Unass. metabol. gene", "Assigned metabol. gene", "other assigned", "other unassigned"))) %>%
  arrange(mapping)

kegg_bar <- kegg_tibble %>%
  mutate(mapping_label = factor(mapping_label, levels = kegg_tibble$mapping_label)) %>%
  ggplot(aes(x = 1, y = gene_count, fill = mapping_label)) +
  geom_col(position = "fill") +
  scale_fill_manual(values = c("#DDA0DD", "#DA70D6", "#555555", "lightgrey")) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin = margin(c(5,5,5,5), unit = "pt"),
  )

#####
# AMG curation


#####
# MORON BAR

#####
# Gene sample prevalence

classification %>%
  filter(contig %in% CDSs_with_metabolism_kegg$contig_id) %>%
  count(Class)

sample_and_country_prevalence <- phold_predictions_with_extensions %>%
  filter(product %in% CDSs_with_metabolism_kegg$product) %>%
  select(cds_id, contig_id, product) %>%
  left_join(., phage_tpm, by = join_by(contig_id == contig)) %>%
  pivot_longer(-c(cds_id, contig_id, product), names_to = "Sample_ID", values_to = "tpm") %>%
  left_join(., metadata[c("Sample_ID", "Country")], by = "Sample_ID") %>%
  filter(tpm > 0) %>% 
  group_by(product) %>%
  summarise(
    number_of_samples = n_distinct(Sample_ID),
    number_of_countries = n_distinct(Country)
    ) %>%
  mutate(
    sample_prevalance = number_of_samples / (ncol(phage_tpm)-1),
    country_prevalence = number_of_countries / n_distinct(metadata$Country, na.rm = TRUE)
  ) %>%
  arrange(desc(sample_prevalance))

sample_and_country_prevalence %>% filter(country_prevalence == 1)

# Of the 5 proteins that occur in all countries, "Levanase" and 
# "Glucosyltransferase" only occur in few samples (11.3% and 13.9%, 
# respectively) and "MazF-like growth inhibitor" and "toxin" are 
# toxin/antitoxins. This leaves only PAPS reductase as protein of particular
# interest.

genes_of_particular_interest <- "PAPS reductase"
genes_of_interest <- c("PAPS reductase", "Glucosyltransferase", "Levanase")

metabolism_tibble <- phold_predictions_with_extensions %>%
  filter(product %in% genes_of_interest) %>%
  count(product, name = "gene_count") %>%
    mutate(asterisk = ifelse(product %in% genes_of_particular_interest, "*", ""),
           product_label = paste0(asterisk, asterisk, product, " (", gene_count, ")", asterisk, asterisk)) %>%
    select(-asterisk) %>%
  arrange(gene_count) %>%
  mutate(product = factor(product, levels = product))

# gene_colors <- c("PAPS reductase" = "#ef8f01",
#                  "Glucosyltransferase" = "#306464",
#                  "Levanase" = "#4A9B9B")
gene_colors <- c("PAPS reductase" = "#ef8f01",
                 "Glucosyltransferase" = "#FFB89A",
                 "Levanase" = "#FF8C69")

goi_bar <- metabolism_tibble %>%
  # mutate(product = factor(product, levels = rev(names(gene_colors)))) %>%
  ggplot(aes(x = 1, y = gene_count, fill = product)) +
  geom_col(position = "fill") +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin = margin(c(5,5,5,5), unit = "pt"),
  ) +
  scale_fill_manual(values = gene_colors,
                    labels = setNames(as.character(metabolism_tibble$product_label),
                                      metabolism_tibble$product)
  )

#####
# GENE PREVALENCE
grene_presence_on_contigs <- phold_predictions_with_extensions %>%
  filter(product %in% metabolism_tibble$product) %>%
  rename(contig = contig_id) %>%
  select(contig, product) %>%
  mutate(present = 1) %>%
  distinct() %>%
  pivot_wider(names_from = product, values_from = present, values_fill = 0)

gene_presence_in_samples <- grene_presence_on_contigs %>%
  pivot_longer(-contig, names_to = "gene", values_to = "present") %>%
  left_join(., phage_tpm, by = "contig") %>%
  mutate(across(-c(contig, gene, present), ~ .x * present)) %>%
  select(-present) %>%
  pivot_longer(-c(contig, gene), names_to = "Sample_ID", values_to = "tpm") %>%
  group_by(gene, Sample_ID) %>%
  mutate(presence = ifelse(sum(tpm) > 0 , 1, 0)) %>%
  ungroup() %>%
  select(-c(contig, tpm)) %>%
  distinct() %>%
  left_join(., metadata[c("Sample_ID", "Country")], by = "Sample_ID") %>%
  group_by(gene) %>%
  ungroup() %>%
  mutate(gene = fct_reorder(gene, presence, .fun = mean, .desc = TRUE))

prevalence_plot <- list()
prevalence_plot$overall <- gene_presence_in_samples %>%
  ggplot(aes(x = reorder(gene, -presence, FUN = mean), y = presence, fill = gene)) +
  geom_bar(stat = "summary", fun = mean) +
  geom_text(
    stat = "summary",
    fun = mean,
    aes(label = paste0(round(after_stat(y*100)), "%")),
    vjust = -0.5,
  ) +  
  labs(x = "Gene", y = "Prevalence") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_fill_manual(values = gene_colors) +
  theme_void() +
  theme(
    plot.margin = margin(r = 5, l = 5, t = 5, b = 5, unit = "pt"),
    axis.text.x = element_text(angle = 90, hjust = 1, margin = margin(t = 5)),
    legend.title = element_text(face = "bold")
  ) 

prevalence_plot$country_facet <- prevalence_plot$overall + 
  facet_wrap(~Country) +
  theme_grey() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1)
  ) 
prevalence_plot$gene_facet <- prevalence_plot$overall + 
  aes(x = Country, y = presence, fill = gene) + facet_wrap(~gene) +
  theme_grey()

#####
# Hosts

host_pie_colors <- c("***Bifidobacterium***" = "#FFDAB9",
                     "***Lactobacillus***" = "#FFA07A",
                     "***Snodgrassella***" = "#FFC300",
                     "***Bombilactobacillus***" = "#ef8f01",
                     "***Gilliamella***" = "#8B4513",
                     "*Frischella*" = "#1C3A3A",
                     "*Commensalibacter*" = "#2A6666",
                     "*Bartonella*" = "#4CB3B3",
                     "*Bombella*" = "#338080",
                     "other" = "#555555",
                     "unknown" = "lightgrey")

hosts_of_genes_tibble <- grene_presence_on_contigs %>%
  pivot_longer(-contig, names_to = "gene") %>%
  filter(value > 0) %>%
  left_join(., classification[c("contig", "Host_group")], by = "contig") %>%
  group_by(gene, Host_group) %>%
  reframe(host_count = n()) %>%
  arrange(gene) %>%
  mutate()

hosts_of_genes_plot_all_genes <- hosts_of_genes_tibble %>%
  mutate(
    # gene = factor(gene, levels = c("PAPS reductase", "Chitinase", "Glucosyltransferase", "Levanase", "PnuC",
    #                                "RimK", "GATase", "RmlC", "NrdD", "CobT", "NMNAT", "Porphyrin synthesis protein")),
    gene = factor(gene, levels = rev(levels(metabolism_tibble$product))),
    Host_group = case_when(Host_group %in% c("Gilliamella", "Lactobacillus", "Bifidobacterium", "Bombilactobacillus", "Snodgrassella") ~ paste0("***", Host_group, "***"),
                           Host_group %in% c("other", "unknown") ~ Host_group,
                           .default = paste0("*", Host_group, "*")),
    Host_group = factor(Host_group, levels = rev(c("unknown", "other", "***Gilliamella***", "***Snodgrassella***", "*Bombella*", "*Bartonella*", "*Frischella*")))
    ) %>%
  ggplot(aes(x = gene, y = host_count, fill = Host_group)) +
  geom_col() +
  scale_fill_manual(values = host_pie_colors) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1),
        legend.text = element_markdown()) +
  labs(x = "Metabolic gene", y = "Host genus count")

hosts_of_genes_plot_goi <- hosts_of_genes_tibble %>%
  filter(gene %in% genes_of_particular_interest) %>%
  mutate(gene = factor(gene, levels = c("PAPS reductase", "Chitinase", "Glucosyltransferase", "Levanase", "PnuC",
                                        "RimK", "GATase", "RmlC", "NrdD", "CobT", "NMNAT", "Porphyrin synthesis protein")),
         Host_group = case_when(Host_group %in% c("Gilliamella", "Lactobacillus", "Bifidobacterium", "Bombilactobacillus", "Snodgrassella") ~ paste0("***", Host_group, "***"),
                                Host_group %in% c("other", "unknown") ~ Host_group,
                                .default = paste0("*", Host_group, "*")),
         Host_group = factor(Host_group, levels = rev(c("unknown", "other", "*Bombella*", "*Bartonella*", "*Frischella*", "***Gilliamella***", "***Snodgrassella***")))) %>%
  rename(`Host genus` = Host_group) %>%
  ggplot(aes(x = gene, y = host_count, fill = `Host genus`)) +
  geom_col() +
  scale_fill_manual(values = host_pie_colors) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=35, vjust=1.1, hjust=1, size = 11),
        plot.margin = margin(c(5,5,5,5), unit = "pt"),
        legend.title = element_text(face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold"),
        legend.text = element_markdown()
  ) +
  labs(x = "Metabolic gene", y = "Genome count")
  
#####
# Disassemble to make a pretty figure

legend_gg <- list()

# PREVALENCE
# Remove geom_text
p <- prevalence_plot$overall
p$layers <- Filter(function(l) !inherits(l$geom, "GeomText"), p$layers)

prev_plot <- p + 
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        plot.margin = margin(c(5,5,5,5), unit = "pt")
  ) +
  coord_flip() +
  scale_y_reverse(expand = expansion(mult = c(0.1, 0))) +
  geom_text(
    stat = "summary",
    fun = mean,
    aes(label = paste0(round(after_stat(y*100)), "%")),
    vjust = 0.5,
    hjust = 1.2,
    size = 5
  )

# PHROG
phrog_with_legend <- phrog_bar_vertical$all +
  guides(
    fill = guide_legend(
      title            = "PHROG category",
      title.theme      = element_text(face = "bold"),
      nrow             = 1,
      byrow            = TRUE,
      reverse          = TRUE,
      direction        = "horizontal"
    )
  ) +
  theme(
    legend.position = "bottom"
  )
legend_gg$phrog <- extract_legend(phrog_with_legend)

# KEGG
kegg_with_legend <- kegg_bar +
  guides(
    fill = guide_legend(
      title            = "KEGG assignment",
      title.theme      = element_text(face = "bold"),
      nrow             = 2,
      byrow            = TRUE,
      reverse          = TRUE,
      direction        = "horizontal"
    )
  ) +
  theme(
    legend.position = "bottom"
  )
legend_gg$kegg <- extract_legend(kegg_with_legend)

# GOI
goi_with_legend <- goi_bar +
  guides(
    fill = guide_legend(
      title            = "Gene product",
      title.theme      = element_text(face = "bold"),
    )
  ) +
  theme(
    legend.position = "right",
    legend.text = element_markdown(),
  )
legend_gg$goi <- extract_legend(goi_with_legend)


# #####
# SAVE FILES
# 
# system("mkdir -p output/R/gene_content")
# system("mkdir -p output/R/genes_pathogens_and_landuse")
# 
# write_delim(phold_predictions_with_extensions, 
#             "output/R/gene_content/phold_predictions_with_extensions_bphage_renamed_genes.tsv",
#             delim = "\t")
# 
# system("mkdir -p output/R/genes_pathogens_and_landuse/phrog_and_kegg")
# 
# ggsave("output/R/genes_pathogens_and_landuse/gene_prevalence.overall.pdf",
#        prevalence_plot$overall, height = 8, width = 8)
# ggsave("output/R/genes_pathogens_and_landuse/gene_prevalence.country_facet.pdf",
#        prevalence_plot$country_facet, height = 10, width = 16)
# ggsave("output/R/genes_pathogens_and_landuse/gene_prevalence.gene_facet.pdf",
#        prevalence_plot$gene_facet, height = 10, width = 18.5)
# 
# for (set in names(phrog_bar_vertical)) {
#   ggsave(paste0("output/R/genes_pathogens_and_landuse/phrog_and_kegg/phrog_bar.vertical.", set, ".pdf"),
#          phrog_bar_vertical[[set]], height = 5.85, width = 0.875)
#   ggsave(paste0("output/R/genes_pathogens_and_landuse/phrog_and_kegg/phrog_bar.horizontal.", set, ".pdf"),
#          phrog_bar_horizontal[[set]], height = 3, width = 6)
# }
# 
# ggsave("output/R/genes_pathogens_and_landuse/phrog_and_kegg/kegg_bar.pdf",
#        kegg_bar, height = 5.85, width = 0.75)
# ggsave("output/R/genes_pathogens_and_landuse/phrog_and_kegg/goi_bar.pdf",
#        goi_bar, height = 5.85, width = 0.75)
# 
# ggsave("output/R/genes_pathogens_and_landuse/phrog_and_kegg/gene_prevalence.overall.nolegend.pdf",
#        prev_plot, height = 5.85, width = 4)
# 
# 
# for (legend in names(legend_gg)) {
#   ggsave(paste0("output/R/genes_pathogens_and_landuse/phrog_and_kegg/legend.", legend, ".pdf"),
#          legend_gg[[legend]], height = 6, width = 6)
# }
# 
# write_delim(hosts_of_genes_tibble, "output/R/genes_pathogens_and_landuse/hosts_of_genes.tsv",
#             delim = "\t")
# 
# write_delim(kegg_and_phold, "output/R/genes_pathogens_and_landuse/kegg_and_phold.tsv",
#             delim = "\t")
# 
# ggsave("output/R/genes_pathogens_and_landuse/hosts_of_genes_all.pdf", hosts_of_genes_plot_all_genes,
#        width = 8, height = 6)
# ggsave("output/R/genes_pathogens_and_landuse/hosts_of_genes_goi.pdf", hosts_of_genes_plot_goi,
#        width = 6, height = 5.2)
# 
# 
# # Update lengths classification table for convenience to avoid backtracking.
# # cobra_refinement_stats <- read.delim("output/core_contig_refinement/cobra_refinement_stats.tsv") %>%
# #   filter(extended_bp > 0)
# # 
# # updated_classification <- classification %>% 
# #   left_join(., cobra_refinement_stats[c("original_name", "refined_length")], by = join_by(contig == original_name)) %>%
# #   mutate(length_kb_after_refinement = case_when(is.na(refined_length) ~ length_kb,
# #                                                 .default = refined_length/1000),
# #          contig_length_refined = case_when(is.na(refined_length) ~ FALSE,
# #                                       .default = TRUE),
# #          .after = length_kb) %>%
# #   select(-refined_length)
# # 
# # saveRDS(updated_classification, "output/R/R_variables/classification.RDS")