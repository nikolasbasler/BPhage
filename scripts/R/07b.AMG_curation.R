library(tidyverse)
library(gggenomes)
metadata <- readRDS("data/metadata.RDS") %>% tibble()

classification <- readRDS("data/classification.RDS")

phold_predictions_with_extensions_bphage_renamed_genes <- read.delim("output/R/gene_content/phold_predictions_with_extensions_bphage_renamed_genes.tsv") %>%
  tibble()

paps_CDSs <- phold_predictions_with_extensions_bphage_renamed_genes %>%
  filter(product == "PAPS reductase") %>%
  select(cds_id) %>%
  unlist(use.names = FALSE)

pois <- c("PAPS reductase", "Chitinase", "Levanase", "PnuC", "Glucosyltransferase")
contigs_with_POI <- list()
for (poi in pois) {
  contigs_with_POI[[poi]] <- phold_predictions_with_extensions_bphage_renamed_genes %>%
    filter(product == poi) %>%
    distinct(contig_id) %>%
    unlist(use.names = FALSE)
}

custom_colors <- c(
  "head and packaging" = "#235050",
  "connector" = "#338080",
  "tail" = "#66CCCC",
  "lysis" = "#8B4513",
  "transcription regulation" = "#ef8f01",
  "integration and excision" = "#FFC300",
  "DNA, RNA and nucleotide metabolism" = "#FFA07A",
  "moron, auxiliary metabolic gene and host takeover" = "#DDA0DD",
  "PAPS reductase" = "#6A0DAD",
  "Chitinase" = "#6A0DAD",
  "Levanase" = "#6A0DAD",
  "PnuC" = "#6A0DAD",     
  "Glucosyltransferase" = "#6A0DAD",
  "unknown function" = "lightgray",
  "empty" = "lightgray",
  "other" = "black",
  "unknown with phrog" = "black"
)

plots_of_genomes <- list()
features <- list()
for (poi in pois) {

  features[[poi]] <- phold_predictions_with_extensions_bphage_renamed_genes %>%
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
    scale_fill_manual(values = custom_colors, breaks = names(custom_colors)) +
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

################################################################################

# For the plots made with the annotations from phage KT1 (see supplements)
# Middelboe, Mathias, Sachia J Traving, Daniel Castillo, Panos G Kalatzis, and 
# Ronnie N Glud. ‘Prophage-Encoded Chitinase Gene Supports Growth of Its 
# Bacterial Host Isolated from Deep-Sea Sediments’. 
# The ISME Journal 19, no. 1 (2025): wraf004. https://doi.org/10.1093/ismejo/wraf004.
# 
# phage_KT1_pharokka <- read_gbk("output/amg_confirmation/chitinase_paper/pharokka_phage_KT1/phage_KT1.gbk") %>%
#   rename(annotation = `function`) %>%
#   mutate(seq_id = "pharokka") %>%
#   # select(seq_id, start, end, product)
#   select(seq_id, start, end, annotation)
# phage_KT1_paper <- read.delim("output/amg_confirmation/chitinase_paper/phage_KT1_annotations_fixed.tsv") %>%
#   tibble() %>%
#   mutate(seq_id = "paper") %>%
#   rename(funct = `Predicted.function`,
#          cons_domain = `Conserved.domain.database.hit.E.value.0.001..accession.`,
#          blast = `BLASTP.most.significant.match..accession..E..value.`) %>%
#   # select(seq_id, start, end, funct, blast, cons_domain)
#   select(seq_id, start, end, cons_domain) %>%
#   rename(annotation = cons_domain)
# 
# feats <- rbind(phage_KT1_pharokka, phage_KT1_paper) %>%
#   mutate(fill_group = ifelse(str_starts(annotation, "Predicted chitinase"), "Chitinase", annotation))
# 
# seqs <- tibble(seq_id = c("pharokka", "paper"),
#                length = c(40573, 40573))
# 
# 
# chitinase_paper_colors <- c(
#   "head and packaging" = "#235050",
#   "connector" = "#338080",
#   "tail" = "#66CCCC",
#   "lysis" = "#8B4513",
#   "transcription regulation" = "#ef8f01",
#   "integration and excision" = "#FFC300",
#   "DNA, RNA and nucleotide metabolism" = "#FFA07A",
#   "moron, auxiliary metabolic gene and host takeover" = "#DDA0DD",
#   "Chitinase" = "#6A0DAD"
# )
# 
# 
# genome_plot <- gggenomes(genes = feats, seqs = seqs) +
#   geom_seq() +
#   geom_seq_label() +
#   geom_gene() +
#   geom_gene(aes(fill=fill_group)) +
#   scale_fill_manual(values = chitinase_paper_colors, breaks = names(chitinase_paper_colors))

# ggsave("output/amg_confirmation/chitinase_paper/phage_KT1_annotations_paper_pharokka.pdf",
#        width = 10, height = 3)

################################################################################

###
# Only keep genes, where
# 1.  They are not the last annotation before one of the contig edges, unless 
#     completeness is 100%
# 2.  Between both sides of the gene and the contig edges, there is at least one
#     other gene assigned to a PHROGS category (unknown genes in between are ok,
#     as long as they are assigned to a PHROG).
# 3.  On one or no side of the gene, the next annotated gene is assigned to a 
#     structural function.

# Only applied to the genomes carrying PAPS reductase, because chitinase fails
# these criteria unanimously and the other GOIs didn't reveal any interesting
# associations downstream.

remove_for_stringency <- list()
remove_for_stringency$`PAPS reductase` <- tibble(
  contig = c(
    "NODE_A10_length_30913_cov_20.162408_NL_19102_spr_mid_d",
    "NODE_A21_length_29443_cov_20.814377_NL_19104_spr_rec_d",
    "NODE_A7_length_26493_cov_31.033502_CH_17692_aut_rec_d",
    "NODE_A8_length_37910_cov_326.687099_PT_19414_sum_rec_d"
    # "NODE_A2_length_44969_cov_51.188675_PT_19409_aut_mid_d" # Don't remove this one! there are 2 fragmented PAPS reductase genes. One is transferred from pharokka, the other one is a pholdseek hit with evalue 1.5e-05
  ),
  reason = c(
    "at_edge_of_incomplete_genome",
    "at_edge_of_incomplete_genome",
    "not_flanked_by_viral_genes",
    "at_edge_of_incomplete_genome"
    # "transferred_pharokka_annotation"
  )
)
remove_for_stringency$Chitinase <- tibble(contig = "", reason = "")
remove_for_stringency$Levanase <- tibble(contig = "", reason = "")
remove_for_stringency$PnuC <- tibble(contig = "", reason = "")
remove_for_stringency$Glucosyltransferase <- tibble(contig = "", reason = "")

paps_completeness_histogram <- classification %>%
  filter(contig %in% contigs_with_POI$`PAPS reductase`,
         !contig %in% remove_for_stringency$`PAPS reductase`$contig) %>%
  ggplot(aes(x = completeness)) +
  geom_histogram(bins = 40) +
  ggtitle("contigs carrying PAPS reductase")

paps_length_histogram <- classification %>%
  filter(contig %in% contigs_with_POI$`PAPS reductase`,
         !contig %in% remove_for_stringency$`PAPS reductase`$contig) %>%
  ggplot(aes(x = length_kb_after_refinement)) +
  geom_histogram(bins = 40) +
  ggtitle("contigs carrying PAPS reductase") +
  labs(x = "contig length (kb)")

goi_presence <- classification %>%
  select(contig)
for (poi in pois) {
  name_of_col <- paste0("GOI_", poi)
  goi_presence <- goi_presence %>%
    mutate(!!name_of_col := case_when(contig %in% contigs_with_POI[[poi]] & !(contig %in% remove_for_stringency[[poi]]$contig) ~ TRUE,
                                      .default = FALSE))
}

genomes_with_goi <- goi_presence %>%
  select(contig, starts_with("GOI_")) %>%
  pivot_longer(-contig) %>%
  filter(value) %>%
  pivot_wider(values_fill = FALSE) %>%
  left_join(., classification, by = "contig")

complete_caudos_with_goi <- classification %>%
  select(-starts_with("GOI_")) %>%
  filter(Class == "Caudoviricetes",
         completeness == 100) %>%
  left_join(., goi_presence, by = "contig")

prop_of_goi_carrying_genomes <- complete_caudos_with_goi %>%
  summarise(across(starts_with("GOI_"), ~ sum(.x, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "GOI", values_to = "abs") %>%
  left_join(
    complete_caudos_with_goi %>%
      summarise(across(starts_with("GOI_"), ~ mean(.x, na.rm = TRUE))) %>%
      pivot_longer(everything(), names_to = "GOI", values_to = "prop"),
    by = "GOI"
  ) %>%
  bind_rows(tibble(GOI = "complete_caudo_genomes", abs = nrow(complete_caudos_with_goi), prop = 1)) %>%
  arrange(desc(GOI))

## TODO: 
# - Mark excluded GOIs in plots?

## Save files

system("mkdir -p output/R/AMG_curation/")

write_delim(remove_for_stringency$`PAPS reductase`, "output/R/AMG_curation/remove_for_stringency.PAPS reductase.tsv",
            delim = "\t")

write_delim(prop_of_goi_carrying_genomes, "output/R/AMG_curation/prop_of_goi_carrying_genomes.tsv",
            delim = "\t")

ggsave("output/R/AMG_curation/histogram_completeness.PAPS reductase.pdf", paps_completeness_histogram,
       width = 4, height = 4)
ggsave("output/R/AMG_curation/histogram_length.PAPS reductase.pdf", paps_length_histogram,
       width = 4, height = 4)

for (poi in pois) {
  hei <- 1 + n_distinct(features[[poi]]$seq_id) / 2
  wid <- max(features[[poi]]$length) / 5000

  ggsave(paste0("output/R/AMG_curation/genome_maps.", poi, ".pdf"),
         plots_of_genomes[[poi]], width = wid, height = hei)
  write_delim(features[[poi]], paste0("output/R/AMG_curation/genome_maps.", poi, ".tsv"),
              delim = "\t")
}


##########################################################################################



# List of AMG-containing contigs for George. Saved to data for convenience

# phold_predictions_with_extensions_bphage_renamed_genes %>%
#   filter(product %in% pois) %>%
#   select(cds_id, product) %>%
#   write_delim("data/AMG_CDSs.tsv", delim = "\t")

# genomes_with_goi %>%
#   select(contig, length_kb_after_refinement, contig_length_refined, lowest_taxon,
#          provirus, gene_count, viral_genes, host_genes, checkv_quality,
#          miuvig_quality, completeness, completeness_method, contamination,
#          Host_genus, Lifestyle_replidec, starts_with("GOI")) %>%
#   write_delim("data/AMG_genomes.tsv", delim = "\t")


# Updating classification table to include presence of GOIs. Only PAPS reductase
# presence is curated.
updated_classification <- classification %>%
  select(!starts_with("GOI_"))
for (poi in names(contigs_with_POI)) {
  name_of_col <- paste0("GOI_", poi)
  updated_classification <- updated_classification %>%
    mutate(!!name_of_col := case_when(contig %in% contigs_with_POI[[poi]] & !(contig %in% remove_for_stringency[[poi]]$contig) ~ TRUE,
                                      .default = FALSE),
           .before = "INPHARED_clustered")

}
# write_delim(updated_classification, "output/R/classification.csv", delim = "\t")


#
