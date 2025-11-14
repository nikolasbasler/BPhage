library(tidyverse)
library(ggpubr)
library(ggtext)
library(patchwork)
library(forcats)
library(gggenomes)

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

vfdb_cds_predictions_extended <- read.delim("output/core_contig_refinement/extended_contigs_phold/sub_db_tophits/vfdb_cds_predictions_long_names.tsv") %>%
  tibble()
vfdb_cds_predictions_unextended <- read.delim("output/annotation/phold_compare_bphage_and_others/sub_db_tophits/vfdb_cds_predictions_long_names.tsv") %>%
  tibble() %>%
  filter(!contig_id %in% vfdb_cds_predictions_extended$contig_id)
vfdb_cds_predictions_with_extensions <- rbind(vfdb_cds_predictions_extended, vfdb_cds_predictions_unextended)

phold_predictions_with_extensions <- phold_annotations_unextended %>%
  filter(str_starts(contig_id, "NODE"),
         !contig_id %in% phold_annotations_extended$contig_id) %>%
  rbind(., phold_annotations_extended) %>%
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

paps_foldseek_results_long_names.tsv <- read.delim("output/annotation/phold_compare_bphage_and_others/paps_foldseek_results_long_names.tsv.gz", header=FALSE) %>%
  tibble() %>%
  rename(
    cds_id = V1,
    hit_protein = V2,
    bitscore = V3,
    fident = V4,
    evalue = V5,
    qStart = V6,
    qEnd = V7,
    qLen = V8,
    tStart = V9,
    tEnd = V10,
    tLen = V11
  )

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

pois <- c("PAPS reductase", "Levanase", "Glucosyltransferase",
          "Chitinase" # This is pre-peer review legacy. Chitinase is not in the list of Vibrant's KOs and based on the surrounding genes, it's most likely not an AMG
          )
contigs_with_POI <- list()
for (poi in pois) {
  contigs_with_POI[[poi]] <- phold_predictions_with_extensions %>%
    filter(product == poi) %>%
    distinct(contig_id) %>%
    unlist(use.names = FALSE)
}

phrog_colors <- c(
  "head and packaging" = "#235050",
  "connector" = "#338080",
  "tail" = "#66CCCC",
  "lysis" = "#8B4513",
  "transcription regulation" = "#ef8f01",
  "integration and excision" = "#FFC300",
  "DNA, RNA and nucleotide metabolism" = "#FFA07A",
  "moron, auxiliary metabolic gene and host takeover" = "#DDA0DD",
  "PAPS reductase" = "#6A0DAD",
  "Levanase" = "#6A0DAD",
  "Glucosyltransferase" = "#6A0DAD",
  "Chitinase" = "#6A0DAD",
  "unknown function" = "lightgray",
  "empty" = "lightgray",
  "other" = "black",
  "unknown with phrog" = "black"
)

plots_of_genomes <- list()
features <- list()
for (poi in pois) {
  
  features[[poi]] <- phold_predictions_with_extensions %>%
    filter(contig_id %in% contigs_with_POI[[poi]]) %>%
    left_join(., classification[c("contig", "length_kb_after_refinement", "completeness", "completeness_method", "virus_score", "fdr")], by = join_by(contig_id == contig)) %>%
    mutate(genome_label = paste0(contig_id, " - CheckV completeness: ", completeness, "%, ", completeness_method, " - geNomad virus score: ", virus_score, " (fdr = ", fdr, ")"),
           length = length_kb_after_refinement * 1000) %>%
    select(contig_id, length, genome_label, start, end, strand, phrog, function., product) %>%
    rename(seq_id = contig_id, funct = function.) %>%
    mutate(
      type = "CDS",
      funct = ifelse(funct == "", "empty", funct),
      fill_group = ifelse(product %in% pois, product, funct),
      fill_group = ifelse(fill_group == "unknown function" & phrog != "No_PHROG", "unknown with phrog", fill_group),
      gene_label = ifelse(product %in% pois, as.character(product), NA_character_),
      contig = seq_id
    )
  
  genome_seqs <- features[[poi]] %>%
    select(seq_id, length, genome_label) %>%
    distinct()
  
  plots_of_genomes[[poi]] <- gggenomes(genes = features[[poi]], seqs = genome_seqs) +
    geom_seq() +
    geom_seq_label(aes(label = genome_label), fontface = "bold") +
    geom_gene(aes(fill=fill_group)) +
    scale_fill_manual(values = phrog_colors, breaks = names(phrog_colors)) +
    geom_gene_tag(aes(label = gene_label, color = "#6A0DAD"),
                  na.rm = TRUE,
                  size = 2.5,
                  angle = 0,
                  nudge_y = 0.14,
                  hjust = 0.5) +
    scale_color_identity(guide = "none") +
    theme(legend.position = "top",
          legend.title = element_text(face = "bold"),
          legend.key.spacing.x = unit(0.5, "cm")
    ) +
    labs(fill = "PHROG category or GOI")
}


paps_phrogs <- c("phrog_2302", "phrog_424", "phrog_33262")
CDSs_with_more_than_one_non_paps_in_top_three_foldseek_hits <- paps_foldseek_results_long_names.tsv %>%
  group_by(cds_id) %>%
  arrange(desc(bitscore)) %>%
  mutate(
    hit_number = row_number(),
    total_hits = n()
  ) %>% 
  ungroup() %>% 
  filter(hit_number <= 3) %>% 
  mutate(hit_phrog = str_extract(hit_protein, "phrog_[0-9]*")) %>% 
  group_by(cds_id, hit_phrog) %>%
  count(hit_phrog) %>%
  filter(
    !hit_phrog %in% paps_phrogs,
    n > 1
    )

###
# Only keep genes, where
# 1.  They are not the last annotation before one of the contig edges, unless 
#     completeness is 100%
# 2.  Between both sides of the gene and the contig edges, there is at least one
#     other gene assigned to a PHROGS category (unknown genes in between are ok,
#     as long as they are assigned to a PHROG).
# 3.  On one or no side of the gene, the next annotated gene is assigned to a 
#     structural function.
# 4.  Max 1 non-paps hit in the top 3 foldseek hits

removed_CDSs <- tibble(
  cds = c(
    "NODE_A10_length_30913_cov_20.162408_NL_19102_spr_mid_d_CDS_0035",
    "NODE_A21_length_29443_cov_20.814377_NL_19104_spr_rec_d_CDS_0001",
    "NODE_A7_length_26493_cov_31.033502_CH_17692_aut_rec_d_CDS_0005",
    "NODE_A8_length_37910_cov_326.687099_PT_19414_sum_rec_d_CDS_0044",
    "NODE_A2_length_44969_cov_51.188675_PT_19409_aut_mid_d_CDS_0039",  # There is another, fragmented PAPS reductase gene on this contig, with a pholdseek hit with evalue 1.5e-05
    "NODE_A3_length_73628_cov_24.519327_BE_16556_sum_rec_d_CDS_0023",
    "NODE_A4_length_61640_cov_572.171532_BE_16556_sum_rec_d_CDS_0091",
    "NODE_A2_length_53539_cov_25.264468_DE_18029_spr_rec_d_CDS_0003",
    "NODE_A3_length_31309_cov_22.437180_DE_18032_sum_mid_d_CDS_0004"
    
    ),
  product = c(
    "PAPS reductase",
    "PAPS reductase",
    "PAPS reductase",
    "PAPS reductase",
    "PAPS reductase",
    "PAPS reductase",
    "PAPS reductase",
    "PAPS reductase",
    "PAPS reductase"
    ),
  reason = c(
    "at_edge_of_incomplete_genome",
    "at_edge_of_incomplete_genome",
    "not_flanked_by_viral_genes",
    "at_edge_of_incomplete_genome",
    "transferred_pharokka_annotation",
    "more_than_1_non_paps_hit_in_top_3_foldseek_hits",
    "more_than_1_non_paps_hit_in_top_3_foldseek_hits",
    "more_than_1_non_paps_hit_in_top_3_foldseek_hits",
    "more_than_1_non_paps_hit_in_top_3_foldseek_hits"
    )
  ) %>%
  rbind(tibble(
    cds = c(
      "NODE_A10_length_37053_cov_15.348794_PT_19407_sum_rec_d_CDS_0041",
      "NODE_A12_length_38214_cov_201.678370_PT_19413_sum_rec_d_CDS_0001",
      "NODE_A1_length_73530_cov_277.772249_RO_17363_spr_mid_d_CDS_0074",
      "NODE_A2_length_32063_cov_15.566779_NL_19104_sum_mid_d_CDS_0029",
      "NODE_A3_length_49905_cov_25.785020_BE_16562_aut_mid_d_CDS_0027",
      "NODE_A3_length_45605_cov_104.639321_PT_19412_sum_mid_d_CDS_0054",
      "NODE_A3_length_35019_cov_23.526272_RO_17363_spr_mid_d_CDS_0057"
      ),
    product = c(
      "Levanase",
      "Levanase",
      "Levanase",
      "Levanase",
      "Levanase",
      "Levanase",
      "Levanase"
      ),
    reason = c(
      "at_edge_of_incomplete_genome",
      "at_edge_of_incomplete_genome",
      "flanked_by_structural_genes",
      "flanked_by_structural_genes",
      "flanked_by_structural_genes",
      "flanked_by_structural_genes",
      "at_edge_of_incomplete_genome"
      )
    )
    ) %>%
  rbind(tibble(
    cds = c(
      "NODE_A12_length_15305_cov_21.290846_NL_19100_sum_mid_d_CDS_0015",
      "NODE_A4_length_43482_cov_12.495703_DE_18029_spr_mid_d_CDS_0009",
      "NODE_A9_length_39063_cov_166.162905_FR_19767_sum_rec_d_CDS_0045"
    ),
    product = c(
      "Glucosyltransferase",
      "Glucosyltransferase",
      "Glucosyltransferase"
      ),
    reason = c(
      "flanked_by_structural_genes",
      "flanked_by_structural_genes",
      "flanked_by_structural_genes"
    )
  )) %>%
  mutate(contig = str_replace(cds, "_CDS_.*", ""), .before = cds)

pois <- c("PAPS reductase", "Levanase", "Glucosyltransferase")
updated_contigs_with_POI <- list()
for (poi in pois) {
  filt_remoded_CDSs <- removed_CDSs %>%
    filter(product == poi)
  updated_contigs_with_POI[[poi]] <- phold_predictions_with_extensions %>%
    filter(product == poi,
           !cds_id %in% filt_remoded_CDSs$cds) %>%
    distinct(contig_id) %>%
    unlist(use.names = FALSE)
}

paps_length_and_completeness_quantiles <- classification %>%
  filter(contig %in% updated_contigs_with_POI$`PAPS reductase`) %>%
  reframe(quantile = names(quantile(completeness)),
          completeness = quantile(completeness),
          length_kb_after_refinement = quantile(length_kb_after_refinement))

paps_completeness_histogram <- classification %>%
  filter(contig %in% updated_contigs_with_POI$`PAPS reductase`) %>%
  ggplot(aes(x = completeness)) +
  geom_histogram(bins = 40) +
  ggtitle("contigs carrying PAPS reductase")

paps_length_histogram <- classification %>%
  filter(contig %in% updated_contigs_with_POI$`PAPS reductase`) %>%
  ggplot(aes(x = length_kb_after_refinement)) +
  geom_histogram(bins = 40) +
  ggtitle("contigs carrying PAPS reductase") +
  labs(x = "contig length (kb)")

goi_presence <- classification %>%
  select(contig) %>%
  mutate(
    `GOI_PAPS reductase` = ifelse(contig %in% updated_contigs_with_POI$`PAPS reductase`, TRUE, FALSE),
    GOI_Levanase = ifelse(contig %in% updated_contigs_with_POI$Levanase, TRUE, FALSE),
    GOI_Glucosyltransferase = ifelse(contig %in% updated_contigs_with_POI$Glucosyltransferase, TRUE, FALSE)
    )

complete_caudos_with_goi <- classification %>%
  select(-starts_with("GOI_")) %>%
  filter(Class == "Caudoviricetes",
         completeness >= 100) %>%
  left_join(., goi_presence, by = "contig")

prop_of_goi_carrying_genomes <- complete_caudos_with_goi %>%
  summarise(across(starts_with("GOI_"), ~ sum(.x, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "genomes", values_to = "count") %>%
  left_join(
    complete_caudos_with_goi %>%
      summarise(across(starts_with("GOI_"), ~ mean(.x, na.rm = TRUE))) %>%
      pivot_longer(everything(), names_to = "genomes", values_to = "proportion"),
    by = "genomes"
  ) %>%
  bind_rows(tibble(genomes = "complete_caudo_genomes", count = nrow(complete_caudos_with_goi), proportion = 1)) %>%
  arrange(desc(genomes))

paps_families <- complete_caudos_with_goi %>%
  filter(`GOI_PAPS reductase`) %>% 
  distinct(Family) %>%
  deframe()

complete_caudo_familiess <- classification %>%
  filter(
    Class == "Caudoviricetes",
    completeness == 100
  ) %>%
  distinct(Family) %>%
  deframe()

paps_family_prevalence <- complete_caudos_with_goi %>%
  filter(Family %in% paps_families) %>%
  group_by(Family) %>%
  summarise(
    contigs_in_family = n(),
    contigs_with_paps = sum(`GOI_PAPS reductase`),
    family_prev = mean(`GOI_PAPS reductase`),
    .groups = "drop") %>%
  mutate(
    contigs_in_family_ignore_ones = ifelse(contigs_in_family == 1, 0, contigs_in_family),
    contigs_with_paps_ignore_ones = ifelse(contigs_in_family == 1, 0, contigs_with_paps)
  ) %>%
  reframe(
    complete_contigs_in_families_with_paps = sum(contigs_in_family),
    complete_contigs_with_paps = sum(contigs_with_paps),
    paps_prevalence_in_families = complete_contigs_with_paps / complete_contigs_in_families_with_paps,
    
    complete_contigs_in_families_with_paps_ignoring_one_counts = sum(contigs_in_family_ignore_ones),
    complete_contigs_with_paps_ignoring_one_counts = sum(contigs_with_paps_ignore_ones),
    paps_prevalence_in_families_ignoring_one_counts = complete_contigs_with_paps_ignoring_one_counts / complete_contigs_in_families_with_paps_ignoring_one_counts
    ) %>%
  pivot_longer(everything(), names_to = "metric") %>%
  rbind(tibble(
    metric = c("families_of_complete_caudos", 
               "families_of_complete_paps_encoding_caudos",
               "families_of_complete_paps_encoding_caudos_prop"),
    value = c(length(complete_caudo_familiess), 
              length(paps_families),
              length(paps_families) / length(complete_caudo_familiess)
              )
  ))

# complete_caudos_with_goi %>%
#   group_by(Family) %>%
#   summarise(family_prev = mean(`GOI_PAPS reductase`),
#             complete_contigs_in_family = n(),
#             weighted_family_prev = family_prev * complete_contigs_in_family) %>% View()
#   reframe(mean_family_prev = mean(family_prev),
#           median_family_prev = median(family_prev))

#####
# MORON BAR

#####
# Gene sample prevalence

classification %>%
  filter(contig %in% CDSs_with_metabolism_kegg$contig_id) %>%
  count(Class)

sample_and_country_prevalence <- phold_predictions_with_extensions %>%
  filter(product %in% CDSs_with_metabolism_kegg$product,
         !cds_id %in% removed_CDSs$cds,
         !product %in% c("toxin", "MazF-like growth inhibitor") # toxin/antitoxin genes
         ) %>%
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

# "PAPS reductase" and "Glucosyltransferase" occur in all countries, but the 
# latter occurs only in few samples (11.1%) and its status as AMG is debatable. 
# This leaves only PAPS reductase as protein of interest.

paps_gene_and_genome_count <- phold_predictions_with_extensions %>%
  filter(product == "PAPS reductase",
         !cds_id %in% removed_CDSs$cds) %>%
  reframe(product = "PAPS reductase", gene_count = n(), genome_count = n_distinct(contig_id))

#####
# PREVALENCE

paps_sample_prev <- sample_and_country_prevalence %>% 
  filter(product == "PAPS reductase") %>% 
  select(sample_prevalance) %>% 
  unlist(use.names = FALSE)

prevalence_plot <- phage_tpm %>%
  pivot_longer(-contig, names_to = "Sample_ID", values_to = "tpm") %>%
  left_join(., metadata[c("Sample_ID", "Country")], by = "Sample_ID") %>%
  mutate(paps_on_contig = ifelse(contig %in% updated_contigs_with_POI$`PAPS reductase` & tpm > 0, TRUE, FALSE)) %>%
  group_by(Sample_ID, Country) %>%
  mutate(paps_in_sample = ifelse(sum(paps_on_contig) > 0, 1, 0)) %>%
  ungroup() %>%
  select(Sample_ID, Country, paps_in_sample) %>%
  distinct() %>%
  ggplot(aes(x = reorder(Country, -paps_in_sample, FUN = mean), y = paps_in_sample)) +
  geom_bar(stat = "summary", fun = mean) +
  geom_text(
    stat = "summary",
    fun = mean,
    aes(label = paste0(round(after_stat(y*100)), "%")),
    vjust = -0.5,
  ) +
  labs(x = "Country", y = "Prevalence") +
  geom_hline(yintercept = paps_sample_prev, linetype = 2) +
  annotate(
    "label",
    x = 7.25, y = 0.7, 
    label = paste0("overall prevalence: \n", round(paps_sample_prev*100), "%"),
    size = 3.5,
    vjust = 0
    ) +
  ggtitle("Sample prevalence of PAPS reductase") +
  theme_minimal()

#####
# Hosts

host_pie_colors <- c("Bifidobacterium" = "#FFDAB9",
                     "Lactobacillus" = "#FFA07A",
                     "Snodgrassella" = "#FFC300",
                     "Bombilactobacillus" = "#ef8f01",
                     "Gilliamella" = "#8B4513",
                     "Frischella" = "#1C3A3A",
                     "Commensalibacter" = "#2A6666",
                     "Bartonella" = "#4CB3B3",
                     "Bombella" = "#338080",
                     "other" = "#555555",
                     "unknown" = "lightgrey")

paps_hosts_tibble <- classification %>%
  filter(contig %in% updated_contigs_with_POI$`PAPS reductase`) %>%
  count(Host_group, name = "host_count") %>%
  mutate(host_label = case_when(
    Host_group %in% c("Gilliamella", "Lactobacillus", "Bifidobacterium", "Bombilactobacillus", "Snodgrassella") ~ paste0("***", Host_group, "*** **(", host_count, ")**"),
    Host_group %in% c("other", "unknown") ~ paste0(Host_group, " (", host_count, ")"),
    .default = paste0("*", Host_group, "* (", host_count, ")")),
    Host_group = factor(Host_group, levels = rev(c("unknown", "other", "Bartonella", "Frischella", "Gilliamella", "Snodgrassella")))
  ) %>%
  arrange(Host_group)

host_group_colors <- tibble(Host_group = names(host_pie_colors), color = host_pie_colors) %>%
  left_join(paps_hosts_tibble, ., by = "Host_group") %>%
  select(host_label, color) %>%
  deframe()

paps_hosts_bar <- paps_hosts_tibble %>%
  mutate(host_label = factor(host_label, levels = paps_hosts_tibble$host_label)) %>%
  ggplot(aes(x = "", y = host_count, fill = host_label)) +
  geom_col() +
  scale_fill_manual(values = host_group_colors) +
  theme_minimal() +
  theme(
    legend.text = element_markdown(),
    legend.title = element_text(face = "bold")
    ) +
  labs(x = NULL, y = "Genome count", fill = "Host group")

#####
# Disassemble to make a pretty figure

legend_gg <- list()

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

# Hosts
host_bar_without_legend <- paps_hosts_bar +
  theme_void() +
  theme(legend.position = "none")
legend_gg$host_left <- extract_legend(paps_hosts_bar)

p <- paps_hosts_bar +
  theme(legend.position = "top") +
  guides(fill = guide_legend(reverse = TRUE))
legend_gg$host_top <- extract_legend(p)

# #####
# SAVE FILES
system("mkdir -p output/R/gene_content/phrog_kegg_and_host")
system("mkdir -p output/R/gene_content/amg_curation")

write_delim(phold_predictions_with_extensions,
            "output/R/gene_content/phold_predictions_with_extensions_bphage_renamed_genes.tsv",
            delim = "\t")

ggsave("output/R/gene_content/gene_prevalence.PAPS reductase.pdf",
       prevalence_plot, height = 4.5, width = 4.5)

for (set in names(phrog_bar_vertical)) {
  ggsave(paste0("output/R/gene_content/phrog_kegg_and_host/phrog_bar.vertical.", set, ".pdf"),
         phrog_bar_vertical[[set]], height = 5.85, width = 0.875)
  ggsave(paste0("output/R/gene_content/phrog_kegg_and_host/phrog_bar.horizontal.", set, ".pdf"),
         phrog_bar_horizontal[[set]], height = 3, width = 6)
}
ggsave("output/R/gene_content/phrog_kegg_and_host/kegg_bar.pdf",
       kegg_bar, height = 5.85, width = 0.75)

for (legend in names(legend_gg)) {
  ggsave(paste0("output/R/gene_content/phrog_kegg_and_host/legend.", legend, ".pdf"),
         legend_gg[[legend]], height = 6, width = 6)
}

write_delim(paps_hosts_tibble, "output/R/gene_content/hosts.PAPS reductase.tsv",
            delim = "\t")

ggsave("output/R/gene_content/phrog_kegg_and_host/hosts.PAPS reductase.pdf",
       host_bar_without_legend, height = 5.85, width = 0.75)
ggsave("output/R/gene_content/hosts.PAPS reductase_with_legend.pdf",
       paps_hosts_bar, height = 5, width = 4)

write_delim(kegg_and_phold, "output/R/gene_content/kegg_and_phold.tsv",
            delim = "\t")

for (gene in names(updated_contigs_with_POI)) {
  write_lines(updated_contigs_with_POI[[gene]], 
              paste0("output/R/gene_content/amg_curation/contigs_with.", gene, ".txt"))
}
write_delim(removed_CDSs, "output/R/gene_content/amg_curation/removed_CDSs.tsv", delim = "\t")
write_delim(prop_of_goi_carrying_genomes, "output/R/gene_content/amg_curation/complete_caudo_genomes_that_carry_goi.tsv", delim = "\t")
write_delim(paps_family_prevalence, "output/R/gene_content/amg_curation/complete_caudo_genomes_paps_family_prevalence.tsv", delim = "\t")

write_delim(paps_length_and_completeness_quantiles, "output/R/gene_content/amg_curation/paps_length_and_completeness_quantiles.tsv", delim = "\t")
ggsave("output/R/gene_content/amg_curation/paps_length_histogram.pdf",
       paps_length_histogram, height = 4, width = 4)
ggsave("output/R/gene_content/amg_curation/paps_completeness_histogram.pdf",
       paps_completeness_histogram, height = 4, width = 4)

for (poi in names(features)) {
  hei <- 1 + n_distinct(features[[poi]]$seq_id) / 2
  wid <- max(features[[poi]]$length) / 5000
  
  ggsave(paste0("output/R/gene_content/amg_curation/genome_maps.", poi, ".pdf"),
         plots_of_genomes[[poi]], width = wid, height = hei)
  write_delim(features[[poi]], paste0("output/R/gene_content/amg_curation/genome_maps.", poi, ".tsv"),
              delim = "\t")
}

# # Update lengths classification table for convenience to avoid backtracking.
# cobra_refinement_stats <- read.delim("output/core_contig_refinement/cobra_refinement_stats.tsv") %>%
#   filter(extended_bp > 0)
# 
# updated_classification <- classification %>%
#   select(-contains("refine")) %>%
#   left_join(., cobra_refinement_stats[c("original_name", "refined_length")], by = join_by(contig == original_name)) %>%
#   mutate(length_kb_after_refinement = case_when(is.na(refined_length) ~ length_kb,
#                                                 .default = refined_length/1000),
#          contig_length_refined = case_when(is.na(refined_length) ~ FALSE,
#                                       .default = TRUE),
#          .after = length_kb) %>%
#   select(-refined_length)

# # Updating classification table to include presence of GOIs. Only PAPS reductase
# # presence is curated.
# updated_classification <- updated_classification %>%
#   select(!starts_with("GOI_"))
# for (poi in names(updated_contigs_with_POI)) {
#   name_of_col <- paste0("GOI_", poi)
#   updated_classification <- updated_classification %>%
#     mutate(!!name_of_col := case_when(contig %in% updated_contigs_with_POI[[poi]] ~ TRUE,
#                                       .default = FALSE),
#            .before = "INPHARED_clustered")
# }

# write_delim(updated_classification, "output/R/classification.csv", delim = "\t")

# List of AMG-containing contigs for George. Saved to data for convenience
# AMG_CDSs <- phold_predictions_with_extensions %>%
#   filter(product %in% pois) %>%
#   select(cds_id, product, contig_id) %>%
#   left_join(., classification[c("contig", "contig_length_refined")], by = join_by(contig_id == contig)) %>%
#   filter(contig_id %in% unlist(updated_contigs_with_POI))
# write_delim(AMG_CDSs, "data/AMG_CDSs.tsv", delim = "\t")
# 
# AMG_genomes <- updated_classification %>%
#   filter(contig %in% AMG_CDSs$contig_id) %>%
#   select(contig, length_kb_after_refinement, contig_length_refined, lowest_taxon,
#          provirus, gene_count, viral_genes, host_genes, checkv_quality,
#          miuvig_quality, completeness, completeness_method, contamination,
#          Host_genus, Lifestyle_replidec, starts_with("GOI"))
# write_delim(AMG_genomes, "data/AMG_genomes.tsv", delim = "\t")

