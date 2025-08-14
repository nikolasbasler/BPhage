library(tidyverse)
library(ggpubr)
library(ggtext)
library(patchwork)
library(forcats)

source("scripts/R/helpers/mixed_helpers.R")

metadata <- readRDS("data/metadata.RDS")
classification <- readRDS("output/R/R_variables/classification.RDS")
present_in_all_countries <- read_lines("data/core_contigs.txt")

phold_predictions_with_extensions <- read.csv("output/R/gene_content/phold_predictions_with_extensions.csv") %>%
  tibble() %>%
  filter(str_starts(contig_id, "NODE")) %>%
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
  left_join(., phold_predictions_with_extensions[c("contig_id", "cds_id", "phrog", "function.", "product")], by = "cds_id")

phage_tpm <- read.csv("output/R/relative_abundance/phage_tpm.csv") %>%
  tibble()

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

kegg_tibble <- phold_predictions_with_extensions %>%
  filter(product %in% genes_with_kegg) %>%
  reframe(placeholder = "placeholder",
          `Assigned metabol. gene` = length(CDSs_with_metabolism_kegg),
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

gene_colors <- rev(c("PAPS reductase" = "#ef8f01", 
                     "Chitinase" = "#8B4513",
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
  mutate(gene = factor(gene, levels = c("PAPS reductase", "Chitinase", "Glucosyltransferase", "Levanase", "PnuC",
                                        "RimK", "GATase", "RmlC", "NrdD", "CobT", "NMNAT", "Porphyrin synthesis protein")),
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
  filter(gene %in% genes_of_interest) %>%
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

write_delim(phold_predictions_with_extensions, 
            "output/R/gene_content/phold_predictions_with_extensions_bphage_renamed_genes.tsv",
            delim = "\t")

system("mkdir -p output/R/genes_pathogens_and_landuse/phrog_and_kegg")

ggsave("output/R/genes_pathogens_and_landuse/gene_prevalence.overall.pdf",
       prevalence_plot$overall, height = 8, width = 8)
ggsave("output/R/genes_pathogens_and_landuse/gene_prevalence.country_facet.pdf",
       prevalence_plot$country_facet, height = 10, width = 16)
ggsave("output/R/genes_pathogens_and_landuse/gene_prevalence.gene_facet.pdf",
       prevalence_plot$gene_facet, height = 10, width = 18.5)

for (set in names(phrog_bar_vertical)) {
  ggsave(paste0("output/R/genes_pathogens_and_landuse/phrog_and_kegg/phrog_bar.vertical.", set, ".pdf"),
         phrog_bar_vertical[[set]], height = 5.85, width = 0.875)
  ggsave(paste0("output/R/genes_pathogens_and_landuse/phrog_and_kegg/phrog_bar.horizontal.", set, ".pdf"),
         phrog_bar_horizontal[[set]], height = 3, width = 6)
}

ggsave("output/R/genes_pathogens_and_landuse/phrog_and_kegg/kegg_bar.pdf",
       kegg_bar, height = 5.85, width = 0.75)
ggsave("output/R/genes_pathogens_and_landuse/phrog_and_kegg/goi_bar.pdf",
       goi_bar, height = 5.85, width = 0.75)

ggsave("output/R/genes_pathogens_and_landuse/phrog_and_kegg/gene_prevalence.overall.nolegend.pdf",
       prev_plot, height = 5.85, width = 4)


for (legend in names(legend_gg)) {
  ggsave(paste0("output/R/genes_pathogens_and_landuse/phrog_and_kegg/legend.", legend, ".pdf"),
         legend_gg[[legend]], height = 6, width = 6)
}

write_delim(hosts_of_genes_tibble, "output/R/genes_pathogens_and_landuse/hosts_of_genes.tsv",
            delim = "\t")

write_delim(kegg_and_phold, "output/R/genes_pathogens_and_landuse/kegg_and_phold.tsv",
            delim = "\t")

ggsave("output/R/genes_pathogens_and_landuse/hosts_of_genes_all.pdf", hosts_of_genes_plot_all_genes,
       width = 8, height = 6)
ggsave("output/R/genes_pathogens_and_landuse/hosts_of_genes_goi.pdf", hosts_of_genes_plot_goi,
       width = 6, height = 5.2)

