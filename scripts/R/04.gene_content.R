library(tidyverse)
library(ggpubr)
library(ggtext)
library(patchwork)
library(forcats)

source("scripts/R/helpers/mixed_helpers.R")

metadata <- readRDS("output/R/R_variables/metadata.RDS")
classification <- readRDS("output/R/R_variables/classification.RDS")
present_in_all_countries <- read_lines("data/core_contigs.txt")

phold_predictions_with_extensions <- read.csv("output/R/gene_content/phold_predictions_with_extensions.csv") %>%
  tibble() %>%
  filter(str_starts(contig_id, "NODE")) %>%
  # These are protein names. Gene names may be more appropriate, but phold seems to refer to "proucts". 
  mutate(
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
    )

kegg_mapping <- read.delim("data/kegg_mapping.tsv", colClasses = "character") %>%
  tibble()

kegg_and_phold <- kegg_mapping %>%
  left_join(., phold_predictions_with_extensions[c("cds_id", "phrog", "function.", "product")], by = "cds_id")

phage_tpm <- read.csv("output/R/relative_abundance/phage_tpm.csv") %>%
  tibble()


#####
# PHROG BARS

moron_bar_colors <- rev(c("lightgrey", "#FFC300", "#555555"))

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
                                   .default = "Other PHROG\ncategories"
                                   )) %>%
  group_by(collapsed_cat) %>%
  summarise(core = sum(count_core),
            all = sum(count_all)) %>%
  mutate(collapsed_cat = factor(collapsed_cat, levels = c("Other PHROG\ncategories", '"moron", AMG\nand host takeover', "unknown function")))

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
    ggplot(aes(x = name, y = value, fill = collapsed_cat)) +
    geom_col(position = "fill") +
    theme_void() +
    theme(
      legend.position = "bottom",
      plot.margin = margin(r = 5, l = 5, t = 5, b = 5),
      legend.title = element_text(face = "bold")
    ) +
    scale_fill_manual(values = moron_bar_colors) +
    labs(fill = "PHROG category") +
    guides(fill = guide_legend(byrow = TRUE, reverse = TRUE))
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
  filter(Pathway_category == "Metabolism" | 
           product %in% c("Chitinase", "GATase", "PnuC")) %>%
  
  filter(!product %in% c("decoy of host sigma70", "MazF-like growth inhibitor",
                         "toxin", "VFDB virulence factor protein")) %>%
  filter(!str_detect(product, "Que")) %>% # This will remove 3 genes. All of them are only present in one sample (the same one for all 3)
  distinct(cds_id) %>% 
  unlist(use.names = FALSE)

genes_with_kegg <- phold_predictions_with_extensions %>% 
  filter(cds_id %in% CDSs_with_metabolism_kegg) %>% 
  distinct(product) %>%
  unlist(use.names = FALSE)

kegg_or_not_kegg <- phold_predictions_with_extensions %>%
  mutate(kegg_mapping = ifelse(product %in% genes_with_kegg, TRUE, FALSE),
         is_moron = ifelse(function. == "moron, auxiliary metabolic gene and host takeover", TRUE, FALSE)) %>%
  reframe(mapped_to_kegg = c("No", "Yes"),
            gene_count = c(sum(is_moron) - sum(kegg_mapping), sum(kegg_mapping))) %>%
  mutate(mapped_to_kegg = paste0(mapped_to_kegg, " (", gene_count, ")"))

kegg_bar <- kegg_or_not_kegg %>%
  mutate(mapped_to_kegg = factor(mapped_to_kegg, levels = kegg_or_not_kegg$mapped_to_kegg)) %>%
  ggplot(aes(x = 1, y = gene_count, fill = mapped_to_kegg)) +
  geom_col(position = "fill") +
  scale_fill_manual(values = c("#555555", "#ef8f01")) +
  theme_void() +
  theme(
    plot.margin = margin(r = 5, l = 5, b = 5, t = 5, unit = "pt"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  ) +
  guides(fill = guide_legend(title = "Mapped to KEGG",
                             label.position = "left",
                             reverse = TRUE))



#####
# MORON BAR

genes_of_interest <- c("Chitinase",
                       "Glucosyltransferase",
                       "Levanase",
                       "PAPS reductase", 
                       "PnuC")

# The 64 sulf genes occur on 60 diffferent genomes. All other genes only appear once on each genome
metabolism_tibble <- phold_predictions_with_extensions %>%
  filter(product %in% genes_with_kegg) %>%
  count(product, name = "gene_count") %>%
  arrange(desc(gene_count)) %>%
  mutate(asterisk = ifelse(product %in% genes_of_interest, "*", ""),
         product_label = paste0(asterisk, asterisk, product, " (", gene_count, ")", asterisk, asterisk)) %>%
  select(-asterisk)

gene_colors <- rev(c("PAPS reductase" = "#8B4513", 
                     "Chitinase" = "#ef8f01",
                     "Glucosyltransferase" = "#FFC300", 
                     "Levanase" = "#FFA07A", 
                     "RimK" = "#66CCCC", 
                     "PnuC" = "#FFDAB9",
                     "GATase" = "#4CB3B3",
                     "RmlC" = "#338080",
                     "NrdD" = "#2A6666",
                     "CobT" = "#235050",
                     "NMNAT" = "#1C3A3A",
                     "Porphyrin synthesis protein" = "#162E2E"))

goi_bar <- metabolism_tibble %>%
  mutate(product_label = factor(product_label, levels = rev(metabolism_tibble$product_label)),
         product = factor(product, levels = rev(metabolism_tibble$product))) %>%
  ggplot(aes(x = 1, y = gene_count, fill = product)) +
  geom_col(position = "fill") +
  theme_void() +
  theme(
    legend.text = element_markdown(),
    plot.margin = margin(r = 5, l = 5, b = 5, t = 5, unit = "pt"),
    legend.title = element_text(face = "bold")
    ) +
  guides(fill = guide_legend(title = "Gene product")) +
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
    aes(label = round(after_stat(y), 2)),
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
# PATCH IT 

# Remove geom_text
p <- prevalence_plot$overall
p$layers <- Filter(function(l) !inherits(l$geom, "GeomText"), p$layers)

prev_plot <- p + 
  theme(legend.position = "none",
        axis.text.x = element_blank(),
  ) +
  coord_flip() +
  geom_text(
    stat = "summary",
    fun = mean,
    aes(label = round(after_stat(y), 2)),
    vjust = 0.5,
    hjust = -0.25
  )

legend_phrog <- extract_legend(phrog_bar_vertical$all)
legend_kegg <- extract_legend(kegg_bar)
# legend_goi <- extract_legend(goi_bar)

pure_phrog <- phrog_bar_vertical$all + theme(legend.position = "none")
pure_kegg <- kegg_bar + theme(legend.position = "none")
# pure_goi <- goi_bar + theme(legend.position = "none")

bar_patch <- pure_phrog + plot_spacer() + pure_kegg + plot_spacer() + goi_bar + plot_spacer() + prev_plot + plot_layout(nrow = 1, widths = c(1.2, 0.5, 1, 0.5, 1, 0.1, 3) )
legends_plot <- legend_phrog / legend_kegg

# #####
# SAVE FILES

system("mkdir -p output/R/genes_pathogens_and_landuse/phrog_and_kegg")

ggsave("output/R/genes_pathogens_and_landuse/gene_prevalence.overall.pdf",
       prevalence_plot$overall, height = 8, width = 8)
ggsave("output/R/genes_pathogens_and_landuse/gene_prevalence.country_facet.pdf",
       prevalence_plot$country_facet, height = 10, width = 16)
ggsave("output/R/genes_pathogens_and_landuse/gene_prevalence.gene_facet.pdf",
       prevalence_plot$gene_facet, height = 10, width = 18.5)

for (set in names(phrog_bar_vertical)) {
  ggsave(paste0("output/R/genes_pathogens_and_landuse/phrog_and_kegg/phrog_bar.vertical.", set, ".pdf"),
         phrog_bar_vertical[[set]], height = 12, width = 6)
  ggsave(paste0("output/R/genes_pathogens_and_landuse/phrog_and_kegg/phrog_bar.horizontal.", set, ".pdf"),
         phrog_bar_horizontal[[set]], height = 3, width = 6)
}

ggsave("output/R/genes_pathogens_and_landuse/phrog_and_kegg/kegg_bar.pdf",
       kegg_bar, height = 12, width = 6)

ggsave("output/R/genes_pathogens_and_landuse/phrog_and_kegg_patch.pdf",
       bar_patch, height = 6, width = 13.5)
ggsave("output/R/genes_pathogens_and_landuse/phrog_and_kegg_legends.pdf",
       legends_plot, height = 3, width = 6)
