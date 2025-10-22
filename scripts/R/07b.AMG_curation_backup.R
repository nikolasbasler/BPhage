library(tidyverse)
library(gggenomes)
metadata <- readRDS("data/metadata.RDS") %>% tibble()

classification <- readRDS("data/classification.RDS")

phold_predictions_with_extensions_bphage_renamed_genes <- read.delim("output/R/gene_content/phold_predictions_with_extensions_bphage_renamed_genes.tsv") %>%
  tibble()

#

# phold_predictions_with_extensions_bphage_renamed_genes %>%
#   filter(product == "PAPS reductase") %>%
#   View()

#

# TODO: put these files to midsave if they turn out to be used!
VIBRANT_annotations <- read.delim("output/amg_confirmation/VIBRANT_bphage_ALL_1kb_all_viruses/VIBRANT_results_bphage_ALL_1kb_all_viruses/VIBRANT_annotations_bphage_ALL_1kb_all_viruses.tsv") %>%
  tibble()
VIBRANT_annotation_coords <- read.delim("output/amg_confirmation/VIBRANT_bphage_ALL_1kb_all_viruses/bphage_ALL_1kb_all_viruses.annotation.coords") %>%
  tibble()
DRAMv_annotations <- read.delim("output/amg_confirmation/dram-v_bphage_AMG_genomes/annotations.tsv") %>%
  mutate(contig = str_replace(X, '__.*$', '')) %>%
  tibble()

phold_v0_foldseek_results_amg_cds <- read.delim("output/amg_confirmation/phold_v0_foldseek_results_amg_cds_long_names.tsv.gz", header=FALSE) %>%
  tibble() %>%
  rename("cds_id" = V1, "hit_protein" = V2, "bitscore" = V3, "fident" = V4, 
         "evalue" = V5, "qStart" = V6, "qEnd" = V7, "qLen" = V8, "tStart" = V9, 
         "tEnd" = V10, "tLen" = V11)
phold_v1_foldseek_results_amg_cds <- read.delim("output/amg_confirmation/phold_v1_foldseek_results_amg_cds_long_names.tsv.gz", header=FALSE) %>%
  tibble() %>%
  rename("cds_id" = V1, "hit_protein" = V2, "bitscore" = V3, "fident" = V4,
         "evalue" = V5, "qStart" = V6, "qEnd" = V7, "qLen" = V8, "tStart" = V9,
         "tEnd" = V10, "tLen" = V11, "alntmscore" = V12, "lddt" = V13)
phold_v1_per_cds_predictions_long_names <- read.delim("output/amg_confirmation/phold_v1_compare_amg_genomes/phold_per_cds_predictions_long_names.tsv") %>%
  tibble()

#

v1 <- phold_predictions_with_extensions_bphage_renamed_genes %>%
  filter(product == "PAPS reductase") %>%
  select(contig_id, cds_id, start, end, strand, product) %>%
  left_join(., phold_v1_foldseek_results_amg_cds, by = "cds_id")


##

paps_CDSs <- phold_predictions_with_extensions_bphage_renamed_genes %>%
  filter(product == "PAPS reductase") %>%
  select(cds_id) %>%
  unlist(use.names = FALSE)

high_confidence_paps <- phold_v1_per_cds_predictions_long_names %>%
  filter(
    cds_id %in% paps_CDSs,
    product == "phosphoadenosine phosphosulfate reductase",
    # annotation_confidence %in% c("high", "medium"),
    annotation_confidence == "high",
    # alntmscore > 0.5, # redundant to annotation_confidence == "high",
    # lddt > 0.8
  )

high_confidence_paps %>% distinct(contig_id) %>% nrow()
high_confidence_paps %>% distinct(cds_id) %>% nrow()

phold_v1_foldseek_results_amg_cds %>%
  filter(cds_id %in% high_confidence_paps$cds_id) %>%
  group_by(cds_id) %>%
  mutate(hit_no = row_number(),
         hit = str_split_i(hit_protein, ":", 1),
         hit_phrog = str_replace(hit, "^.*_phrog", "phrog")
         ) %>%
  ungroup() %>%
  filter(hit_no <= 10) %>%
  count(hit_phrog) %>%
  arrange(desc(n))

# phold_v0_foldseek_results_amg_cds %>%
#   filter(cds_id %in% paps_CDSs) %>%
#   group_by(cds_id) %>%
#   mutate(hit_no = row_number(),
#          hit_phrog = str_extract(hit_protein, "phrog_[0-9]*")
#   ) %>%
#   ungroup() %>%
#   filter(hit_no <= 10) %>%
#   count(hit_phrog) %>%
#   arrange(desc(n))

high_confidence_paps_where_top_5_are_paps <- phold_v1_foldseek_results_amg_cds %>%
  filter(cds_id %in% high_confidence_paps$cds_id) %>%
  group_by(cds_id) %>%
  mutate(hit_no = row_number(),
         hit_phrog = str_split_i(hit_protein, ":", 1),
         ) %>%
  filter(hit_no <= 10) %>%
  mutate(total_top_hits = n()) %>%
  filter(str_detect(hit_phrog, "phrog_424|phrog_2302|phrog_12882|phrog_33262")) %>%
  mutate(hits_left = n()) %>%
  ungroup() %>%
  # filter(hits_left == total_top_hits) %>%
  filter(hits_left == 10) %>%
  reframe(cds_id = unique(cds_id)) %>%
  mutate(contig = str_replace(cds_id, "_CDS_.*", ""))

# high_confidence_paps_where_top_5_are_paps <- phold_v0_foldseek_results_amg_cds %>%
#   group_by(cds_id) %>%
#   mutate(hit_no = row_number(),
#          hit_phrog = str_split_i(hit_protein, ":", 1),
#   ) %>%
#   filter(hit_no <= 10) %>%
#   mutate(total_top_hits = n()) %>%
#   filter(str_detect(hit_phrog, "phrog_424|phrog_2302|phrog_12882|phrog_33262")) %>%
#   mutate(hits_left = n()) %>%
#   ungroup() %>%
#   # filter(hits_left == total_top_hits) %>%
#   filter(hits_left == 10) %>%
#   reframe(cds_id = unique(cds_id)) %>%
#   mutate(contig = str_replace(cds_id, "_CDS_.*", ""))

# high_confidence_paps_where_top_5_are_paps %>% distinct(contig) %>% nrow()

#

classification %>%
  filter(if_any(starts_with("GOI"), ~ . == TRUE)) %>%
  count(Class)

genomes_with_goi <- classification %>%
  filter(if_any(starts_with("GOI"), ~ . == TRUE))

#

completeness_histogram <- list()
length_histogram <- list()
viral_and_host_genes <- list()
gois <- c("GOI_cysH", "GOI_chitinase", "GOI_levanase", "GOI_pnuc", "GOI_glucosyltransferase")
for (goi in gois) {
  
  filt_classi <- classification %>%
    filter(.data[[goi]])
  
  genomes_with_goi <- nrow(filt_classi)
  
  viral_and_host_genes[[goi]] <- filt_classi %>%
    reframe(
      genomes = genomes_with_goi,
      genes = sum(gene_count),
      viral_genes = sum(viral_genes),
      host_genes = sum(host_genes)
      )
  
  
  completeness_histogram[[goi]] <- filt_classi %>%
    ggplot(aes(x = completeness)) +
    geom_histogram(bins = genomes_with_goi) +
    ggtitle(goi)
  
  length_histogram[[goi]] <- filt_classi %>%
    ggplot(aes(x = length_kb)) +
    geom_histogram(bins = genomes_with_goi) +
    ggtitle(goi)
  
}
  
# completeness_histogram$GOI_cysH
# length_histogram$GOI_cysH
# viral_and_host_genes$GOI_cysH


pois <- c("PAPS reductase", "Chitinase", "Levanase", "PnuC", "Glucosyltransferase")

distances_to_contig_edges <- phold_predictions_with_extensions_bphage_renamed_genes %>%
  filter(product %in% pois) %>%
  left_join(., classification[c("contig", "length_kb_after_refinement", "completeness", "completeness_method")], by = join_by(contig_id == contig)) %>%
  mutate(contig_length = length_kb_after_refinement * 1000) %>%
  select(contig_id, cds_id, product, start, end, contig_length, completeness, completeness_method) %>%
  rowwise() %>%
  mutate(distance_of_goi_to_contig_edge = min(contig_length - max(start, end), min(start, end))) %>%
  select(cds_id, distance_of_goi_to_contig_edge)

contigs_with_POI <- list()
for (poi in pois) {
  contigs_with_POI[[poi]] <- phold_predictions_with_extensions_bphage_renamed_genes %>%
    filter(product == poi) %>%
    distinct(contig_id) %>%
    unlist(use.names = FALSE)
}

###

classification %>%
  filter(Class == "Caudoviricetes",
         completeness == 100) %>%
  summarise(complete_caudos = n(),
            with_paps = sum(GOI_cysH)) %>%
  mutate(proportion  = with_paps / complete_caudos)


###

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
  "VIBRANT_AMG" = "red",
  "DRAMv_aux" = "red",
  "unknown function" = "lightgray",
  "empty" = "lightgray",
  "other" = "lightgray",
  "VIBRANT_unknown" = "lightgray",
  "DRAMv_other" = "lightgray",
  "VIBRANT_other" = "black"
)


## TEST THIS!
vibrant_features <- left_join(VIBRANT_annotations, VIBRANT_annotation_coords, by = join_by(protein == CDS)) %>%
  left_join(., classification[c("contig", "length_kb_after_refinement")], by = join_by(scaffold == contig)) %>%
  
  # rowwise() %>%
  # mutate(
  #   v_score = {
  #     vals <- c_across(contains("v.score"))         # collect the row's values
  #     if (all(is.na(vals)))                       # if all NA -> keep NA
  #       NA
  #     else
  #       as.character(max(vals, na.rm = TRUE))                   # otherwise take max ignoring NAs
  #   }
  # ) %>%
  # ungroup() %>%
  
  mutate(length = length_kb_after_refinement * 1000,
         genome_label = "VIBRANT",
         strand = case_when(strand == 1 ~ "+", 
                            strand == -1 ~ "-",
                            .default = NA),
         fill_group = case_when(AMG == "AMG" ~ "VIBRANT_AMG", 
                                funct == "hypothetical protein" ~ "VIBRANT_unknown",
                                .default = "VIBRANT_other"),
         seq_id = paste0(scaffold, "_VIBRANT"),
         type = "CDS",
         product = NA,
         phrog = NA,
         gene_label = NA,
         # gene_label = v_score,
         start_correct = ifelse(strand == "-", end, start),
         end_correct = ifelse(strand == "-", start, end),
         distance_of_goi_to_contig_edge = NA
         ) %>% 
  mutate(start = start_correct,
         end = end_correct) %>%
  select(seq_id, length, genome_label, start, end, strand, phrog, funct, product, type, fill_group, gene_label, scaffold, distance_of_goi_to_contig_edge)


# dram_features <- DRAMv_annotations %>%
#   mutate(scaffold = str_replace(X, '__.*$', '')) %>%
#   left_join(., classification[c("contig", "length_kb_after_refinement")], by = join_by(scaffold == contig)) %>%
#   mutate(length = length_kb_after_refinement * 1000,
#          genome_label = "DRAMv",
#          strand = case_when(strandedness == 1 ~ "+", 
#                             strandedness == -1 ~ "-",
#                             .default = NA),
#          fill_group = ifelse(auxiliary_score < 4, "DRAMv_aux", "DRAMv_other"),
#          seq_id = paste0(scaffold, "_DRAMv"),
#          type = "CDS",
#          product = NA,
#          phrog = NA,
#          gene_label = NA,
#          start = ifelse(strand == "-", end_position, start_position),
#          end = ifelse(strand == "-", start_position, end_position),
#          distance_of_goi_to_contig_edge = NA,
#          funct = NA
#          ) %>%
#   select(seq_id, length, genome_label, start, end, strand, phrog, funct, product, type, fill_group, gene_label, scaffold, distance_of_goi_to_contig_edge)

plots_of_genomes <- list()
features <- list()
features_with_vibrant <- list()
features_with_dram <- list()
for (poi in pois) {

  features[[poi]] <- phold_predictions_with_extensions_bphage_renamed_genes %>%
    filter(contig_id %in% contigs_with_POI[[poi]]) %>%
    left_join(., classification[c("contig", "length_kb_after_refinement", "completeness", "completeness_method", "virus_score", "fdr")], by = join_by(contig_id == contig)) %>%
    left_join(., distances_to_contig_edges, by = "cds_id") %>%
    mutate(genome_label = paste0(contig_id, " - CheckV completeness: ", completeness, "%, ", completeness_method, " - geNomad virus score: ", virus_score, " (fdr = ", fdr, ")"),
           length = length_kb_after_refinement * 1000) %>%
    select(contig_id, length, genome_label, start, end, strand, phrog, function., product, distance_of_goi_to_contig_edge) %>%
    rename(seq_id = contig_id, funct = function.) %>%
    mutate(
      type = "CDS",
      funct = ifelse(funct == "", "empty", funct),
      fill_group = ifelse(product %in% pois, product, funct),
      fill_group = factor(fill_group, levels = names(custom_colors)),
      gene_label = ifelse(product %in% pois, as.character(product), NA_character_),
      contig = seq_id
      )
  
  features_with_vibrant[[poi]] <- vibrant_features %>%
    rename(contig = scaffold) %>%
    filter(contig %in% features[[poi]]$seq_id) %>%
    rbind(features[[poi]], .) %>%
    arrange(seq_id) %>%
    group_by(contig, start) %>%
    mutate(has_start_with_goi = ifelse(length(na.omit(gene_label)) == 0, FALSE, TRUE)) %>%
    group_by(contig, end) %>%
    mutate(has_end_with_goi = ifelse(length(na.omit(gene_label)) == 0, FALSE, TRUE)) %>%
    ungroup() %>%
    mutate(
      # gene_label = case_when(genome_label == "VIBRANT" & (has_start_with_goi | has_end_with_goi) ~ funct,
      #                        .default = gene_label),
      gene_label = ifelse(genome_label == "VIBRANT" & (has_start_with_goi | has_end_with_goi), funct, gene_label),
      gene_label_color = ifelse(genome_label == "VIBRANT", "#444444", "#6A0DAD")
      )
  # 
  # features_with_dram[[poi]] <- dram_features %>%
  #   rename(contig = scaffold) %>%
  #   filter(contig %in% features[[poi]]$seq_id) %>%
  #   rbind(features[[poi]], .) %>%
  #   arrange(seq_id) %>%
  #   group_by(contig, start) %>%
  #   mutate(has_start_with_goi = ifelse(length(na.omit(gene_label)) == 0, FALSE, TRUE)) %>%
  #   group_by(contig, end) %>%
  #   mutate(has_end_with_goi = ifelse(length(na.omit(gene_label)) == 0, FALSE, TRUE)) %>%
  #   ungroup() %>%
  #   mutate(gene_label = ifelse(genome_label == "DRAMv" & (has_start_with_goi | has_end_with_goi), funct, gene_label),
  #          gene_label_color = ifelse(genome_label == "DRAMv", "#444444", "#6A0DAD"))
  # 
  genome_seqs <- features[[poi]] %>%
    select(seq_id, length, genome_label) %>%
    distinct()
  genome_seqs_with_vibrant <- features_with_vibrant[[poi]] %>%
    select(seq_id, length, genome_label) %>%
    distinct()
  # genome_seqs_with_dram <- features_with_dram[[poi]] %>%
  #   select(seq_id, length, genome_label) %>%
  #   distinct()
  
  links_phold_vibrant <- features_with_vibrant[[poi]] %>%
    filter(has_start_with_goi | has_end_with_goi) %>%
    select(seq_id, start, end, contig) %>%
    mutate(is_vibrant = ifelse(str_detect(seq_id, "VIBRANT"), TRUE, FALSE)) %>%
    group_by(contig, is_vibrant) %>%
    mutate(idx = row_number()) %>%
    ungroup() %>%
    pivot_wider(id_cols = c(idx, contig), names_from = is_vibrant, values_from = c(seq_id, start, end)) %>%
    rename(seq_id = seq_id_FALSE, seq_id2 = seq_id_TRUE,
           start = start_FALSE, start2 = start_TRUE,
           end = end_FALSE, end2 = end_TRUE) %>%
    select(-c(idx, contig))
  
  # links_phold_dram <- features_with_dram[[poi]] %>%
  #   filter(has_start_with_goi | has_end_with_goi) %>%
  #   select(seq_id, start, end, contig) %>%
  #   mutate(is_dram = ifelse(str_detect(seq_id, "DRAM"), TRUE, FALSE)) %>%
  #   group_by(contig, is_dram) %>%
  #   mutate(idx = row_number()) %>%
  #   ungroup() %>% 
  #   pivot_wider(id_cols = c(idx, contig), names_from = is_dram, values_from = c(seq_id, start, end)) %>%
  #   rename(seq_id = seq_id_FALSE, seq_id2 = seq_id_TRUE,
  #          start = start_FALSE, start2 = start_TRUE,
  #          end = end_FALSE, end2 = end_TRUE) %>%
  #   select(-c(idx, contig))

  # plots_of_genomes[[poi]] <- gggenomes(genes = features[[poi]], seqs = genome_seqs) +
  plots_of_genomes[[poi]] <- gggenomes(genes = features_with_vibrant[[poi]], seqs = genome_seqs_with_vibrant, links = links_phold_vibrant) +
  # plots_of_genomes[[poi]] <- gggenomes(genes = features_with_dram[[poi]], seqs = genome_seqs_with_dram, links = links_phold_dram) +
    geom_link() +
    geom_seq() +
    geom_seq_label(aes(label = genome_label), fontface = "bold") +
    geom_gene(aes(fill=fill_group)) +
    scale_fill_manual(values = custom_colors, breaks = names(custom_colors)) +
    geom_gene_tag(aes(label = gene_label, color = gene_label_color),
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
  plots_of_genomes[[poi]]

}


paps_coords <- features$`PAPS reductase` %>%
  filter(product == "PAPS reductase") %>%
  select(seq_id, start, end)

DRAMv_annotations %>%
  mutate(seq_id = contig,
         start = ifelse(strandedness == -1, end_position, start_position),
         end = ifelse(strandedness == -1, start_position, end_position)
         ) %>% 
  left_join(paps_coords, ., by = c("seq_id", "start", "end")) 

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
# 
# ggsave("output/amg_confirmation/chitinase_paper/phage_KT1_annotations_paper_pharokka.pdf",
#        width = 10, height = 3)

################################################################################

###
## CONTINUE HERE. DO THIS FOR THE OTHER GENES AS WELL (OR NOT, BECAUSE NOTHING
## INTERESTING IS FOUND WITH THE OTHER ONES ANYWAY)

# Only keep genes, where
# 1.  They are not the last annotation before one of the contig edges, unless 
#     completeness is 100%
# 2.  Between both sides of the gene and the contig edges, there is at least one
#     other gene assigned to a PHROGS category (unkown genes in between are ok).
# 3.  On one or no side of the gene, the next annotated gene is a structural
#     gene.

phold_predictions_with_extensions_bphage_renamed_genes %>%
  filter(product == "PAPS reductase") %>% View()
  summarise(median(evalue, na.rm = TRUE))

phold_predictions_with_extensions_bphage_renamed_genes %>%
  filter(product == "PAPS reductase") %>% 
  ggplot(aes(x = evalue)) +
  geom_histogram()

phold_predictions_with_extensions_bphage_renamed_genes %>%
  filter(product == "PAPS reductase",
         evalue > 1e-6) %>%
  select(contig_id)

remove_for_stringency <- list()
remove_for_stringency$`PAPS reductase` <- tibble(
  contig = c(
    "NODE_A10_length_30913_cov_20.162408_NL_19102_spr_mid_d",
    "NODE_A21_length_29443_cov_20.814377_NL_19104_spr_rec_d",
    "NODE_A25_length_32356_cov_23.103690_NL_19104_aut_rec_d",
    "NODE_A2_length_53539_cov_25.264468_DE_18029_spr_rec_d",
    "NODE_A4_length_36641_cov_59.766710_BE_16557_sum_mid_d",
    "NODE_A7_length_26493_cov_31.033502_CH_17692_aut_rec_d",
    "NODE_A8_length_39429_cov_27.847987_BE_16562_aut_mid_d",
    "NODE_A8_length_37910_cov_326.687099_PT_19414_sum_rec_d",
    
    "NODE_A2_length_44969_cov_51.188675_PT_19409_aut_mid_d",
    "NODE_A4_length_44865_cov_65.538894_PT_19412_spr_rec_d",
    "NODE_A2_length_44969_cov_51.188675_PT_19409_aut_mid_d",
    "NODE_A8_length_37910_cov_326.687099_PT_19414_sum_rec_d"
  ),
  reason = c(
    "at_edge_of_incomplete_genome",
    "at_edge_of_incomplete_genome",
    "not_flanked_by_viral_genes",
    "not_flanked_by_viral_genes",
    "not_flanked_by_viral_genes",
    "not_flanked_by_viral_genes",
    "not_flanked_by_viral_genes",
    "at_edge_of_incomplete_genome",
    "eval_above_1e-6",
    "eval_above_1e-6",
    "transferred_pharokka_annotation",
    "transferred_pharokka_annotation"
  )
)
remove_for_stringency$Chitinase <- tibble()
remove_for_stringency$Levanase <- tibble()
remove_for_stringency$PnuC <- tibble()
remove_for_stringency$Glucosyltransferase <- tibble()

high_confidence_paps_where_top_5_are_paps %>%
  filter(!cds_id %in% remove_for_stringency$`PAPS reductase`$contig) %>%
  distinct(contig)

is_vibrant_amg <- list()
is_vibrant_amg$`PAPS reductase` <- features_with_vibrant$`PAPS reductase` %>%
  filter(!is.na(gene_label),
         fill_group == "VIBRANT_AMG") %>%
  select(contig, start, end, strand, gene_label) %>%
  distinct()

is_vibrant_amg$Chitinase <- tibble()
is_vibrant_amg$Levanase <- tibble()
is_vibrant_amg$PnuC <- tibble()
is_vibrant_amg$Glucosyltransferase <- tibble()

## Save files

system("mkdir -p output/R/AMG_curation/gene_context/")
write_delim(high_confidence_paps_where_top_5_are_paps, "output/R/AMG_curation/gene_context/high_confidence_paps.tsv",
            delim = "\t")

for (poi in pois) {

#   # Without VIBRANT
#   # hei <- 1 + n_distinct(features[[poi]]$seq_id) / 2
#   # wid <- max(features[[poi]]$length) / 5000
#   # 
#   # ggsave(paste0("output/R/AMG_curation/gene_context/genome_maps.", poi, ".pdf"),
#   #        plots_of_genomes[[poi]], width = wid, height = hei)
#   # write_delim(features[[poi]], paste0("output/R/AMG_curation/gene_context/genome_maps.", poi, ".tsv"),
#   #             delim = "\t")
#   
#   # With VIBRANT
#   hei <- 1 + n_distinct(features_with_vibrant[[poi]]$seq_id) / 2
#   wid <- max(features_with_vibrant[[poi]]$length) / 5000
#   ggsave(paste0("output/R/AMG_curation/gene_context/genome_maps.with_vibrant.", poi, ".pdf"),
#          plots_of_genomes[[poi]], width = wid, height = hei, limitsize = FALSE)
#   write_delim(features_with_vibrant[[poi]], paste0("output/R/AMG_curation/gene_context/genome_maps.with_vibrant.", poi, ".tsv"),
#               delim = "\t")
#   
#   ##
#   
  write_delim(remove_for_stringency[[poi]], paste0("output/R/AMG_curation/gene_context/remove_for_stringency.", poi, ".tsv"),
              delim = "\t")
#   write_delim(is_vibrant_amg[[poi]], paste0("output/R/AMG_curation/gene_context/is_vibrant_amg.", poi, ".tsv"),
#               delim = "\t")
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
# 



#
