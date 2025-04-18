library(lme4)
library(lmerTest)
library(DHARMa)
library(ggrepel)
library(forcats)
library(tidyverse)
library(patchwork)
library(gMCPLite)

source("scripts/R/helpers/mixed_helpers.R")

metadata <- readRDS("output/R/R_variables/metadata.RDS") %>%
  mutate(Hive_ID = as.character(Hive_ID))

cropland_fraction <- read.csv("data/land_cover_results.csv") %>%
  tibble() %>%
  mutate(cropland_fraction = cropland_fraction / 100) %>%
  rename(cropland_fraction_2k_radius = cropland_fraction) %>%
  arrange(cropland_fraction_2k_radius)
# make FAOSTAT_added_data:
source("scripts/R/helpers/FAOstat_table.R")

cropland_and_FAO <- FAOSTAT_added_data %>%
  select(Country, Item, est_use_in_2k_radius) %>%
  pivot_wider(id_cols = Country, values_from = est_use_in_2k_radius, names_from = Item) %>%
  # left_join(cropland_fraction[c("Country", "cropland_fraction_2k_radius")], ., by = "Country") %>%
  left_join(FAOSTAT_added_data[c("Country", "cropland_fraction_2k_radius", "Cropland_in_2km_radius")],. , by = "Country") %>%
  distinct() %>%
  select_if(~ !any(is.na(.)))


pathogen_data <- read_excel("data/GlobalBGOOD_WP1_Tier1_Scien.xlsx", skip = 1) %>%
  rename(BGOOD_sample_code = Sample_ID) %>%
  mutate(CBPV = ifelse(str_detect(CBPV, ">40,00"), "41", CBPV),
         BQCV = as.character(BQCV))


pathogens_Cts <- c("DWV A", "DWV B", "ABPV", "CBPV", "BQCV", "SBV", "EFB",
                   "AFB", "N. apis", "N. ceranae", "DWV A")

pathogen_ct <- pathogen_data %>% 
  filter(Years == 2020) %>%
  select(BGOOD_sample_code, all_of(pathogens_Cts)) %>%
  left_join(metadata[c("BGOOD_sample_code", "Bee_pool")], ., by = "BGOOD_sample_code") %>%
  distinct() %>%
  select(-BGOOD_sample_code) %>%
  pivot_longer(-Bee_pool, names_to = "pathogen", values_to = "Ct") %>%
  filter(!is.na(Bee_pool)) %>% 
  mutate(Ct = str_replace_all(Ct, "negative", "41")) %>%
  mutate(Ct = as.numeric(Ct)) %>% 
  filter(!is.na(Ct))


phage_tpm <- read.csv("output/R/relative_abundance/phage_tpm.csv") %>%
  tibble()

phold_predictions_with_extensions <- read.csv("output/R/gene_content/phold_predictions_with_extensions.csv") %>%
  tibble() %>%
  filter(str_starts(contig_id, "NODE"))

kegg_mapping <- read.delim("data/kegg_mapping.tsv", colClasses = "character") %>%
  tibble()

kegg_and_phold <- kegg_mapping %>%
  left_join(., phold_predictions_with_extensions[c("cds_id", "phrog", "function.", "product")], by = "cds_id")

CDSs_with_metabolism_kegg <- kegg_and_phold %>%
  filter(Pathway_category == "Metabolism" |
           product %in% c("chitinase", "glutamine amidotransferase",
                          "PnuC-like nicotinamide mononucleotide transport")) %>%
  filter(!product %in% c("decoy of host sigma70", "MazF-like growth inhibitor",
                         "toxin", "VFDB virulence factor protein")) %>%
  filter(!str_detect(product, "Que")) %>% # This will remove 3 genes. All of them are only present in one sample (the same one for all 3)
  distinct(cds_id) %>%
  unlist(use.names = FALSE)

genes_with_kegg <- phold_predictions_with_extensions %>%
  filter(cds_id %in% CDSs_with_metabolism_kegg) %>%
  distinct(product) %>%
  unlist(use.names = FALSE)

grene_presence_on_contigs <- phold_predictions_with_extensions %>%
  filter(product %in% genes_with_kegg) %>%
  rename(contig = contig_id) %>%
  select(contig, product) %>%
  mutate(present = 1) %>%
  distinct() %>%
  pivot_wider(names_from = product, values_from = present, values_fill = 0)

gene_tpm <- grene_presence_on_contigs %>%
  pivot_longer(-contig, names_to = "gene", values_to = "present") %>%
  left_join(., phage_tpm, by = "contig") %>%
  mutate(across(-c(contig, gene, present), ~ .x * present)) %>%
  select(-present) %>%
  pivot_longer(-c(contig, gene), names_to = "Sample_ID", values_to = "tpm") %>%
  group_by(gene, Sample_ID) %>%
  mutate(tpm = sum(tpm)) %>%
  ungroup() %>%
  select(-contig) %>%
  distinct()

genes_of_interest <- c("chitinase",
                       "glucosyltransferase",
                       "levanase",
                       "phosphoadenosine phosphosulfate reductase", 
                       "PnuC-like nicotinamide mononucleotide transport")

coeffs <- list()

pathogens_of_interest <- c("DWV B", "BQCV", "SBV")
#####
# GENE PRESENCE VS CT
test_tibble_gene_presence <- list()
model_gene_presence <- list()
coeffs$gene_presence <- tibble()
for (goi in genes_of_interest) {
  for (poi in pathogens_of_interest) {
    test_tibble_gene_presence[[goi]][[poi]] <- gene_tpm %>%
      left_join(., metadata[c("Sample_ID", "Bee_pool", "Country", "Hive_ID", "Season", "Gut_part")], by = "Sample_ID") %>%
      left_join(., pathogen_ct, by = "Bee_pool", relationship = "many-to-many") %>%
      filter(gene == goi,
             pathogen == poi) %>%
      mutate(gene_presence = ifelse(tpm > 0, 1, 0))

    model_gene_presence[[goi]][[poi]] <- glmer(
      gene_presence ~ Ct + Season + Gut_part + ( 1 | Hive_ID ),
      data = test_tibble_gene_presence[[goi]][[poi]],
      family = binomial)
    
    has_convergence_issues <- FALSE
    messages <- model_gene_presence[[goi]][[poi]]@optinfo$conv$lme4$messages
    if (!is.null(messages) & any(str_detect(messages, "failed to converge"))) {
      has_convergence_issues <- TRUE
    }
    if (has_convergence_issues) {
      print(paste0(goi," - ", poi, " - Model didn't converge. Removed from the list."))
      model_gene_presence[[goi]][[poi]] <- NULL
    } else {
      coeffs$gene_presence <- summary(model_gene_presence[[goi]][[poi]])$coefficients %>%
        as.data.frame() %>%
        rownames_to_column("metric") %>%
        tibble() %>%
        mutate(pathogen = poi, .before = metric) %>%
        mutate(gene = goi, .before = metric) %>%
        mutate(Item = metric, .before = metric) %>%
        mutate(singular = ifelse(isSingular(model_gene_presence[[goi]][[poi]]), TRUE, FALSE)) %>%
        rbind(coeffs$gene_presence)
    }
  }
}

#####
# LOG TPM VS CT
test_tibble_abundances <- list()
model_abundances <- list()
coeffs$abundances <- tibble()
for (goi in genes_of_interest) {
  for (poi in pathogens_of_interest) {
    test_tibble_abundances[[goi]][[poi]] <- gene_tpm %>%
      left_join(., metadata[c("Sample_ID", "Bee_pool", "Country", "Hive_ID", "Season", "Gut_part")], by = "Sample_ID") %>%
      left_join(., pathogen_ct, by = "Bee_pool", relationship = "many-to-many") %>%
      filter(gene == goi,
             pathogen == poi) %>%
      mutate(log_tpm = log10(tpm)) %>%
      filter(!is.infinite(log_tpm))
    

    model_abundances[[goi]][[poi]] <- lmer(Ct ~ log_tpm + Gut_part + Season +
                                             (1 | Hive_ID ), 
                                           data = test_tibble_abundances[[goi]][[poi]])
    
    has_convergence_issues <- FALSE
    messages <- model_abundances[[goi]][[poi]]@optinfo$conv$lme4$messages
    if (!is.null(messages) & any(str_detect(messages, "failed to converge"))) {
      has_convergence_issues <- TRUE
    }
    if (has_convergence_issues) {
      print(paste0(goi," - ", poi, " - Model didn't converge. Removed from the list."))
      model_abundances[[goi]][[poi]] <- NULL
    } else {
      coeffs$abundances <- summary(model_abundances[[goi]][[poi]])$coefficients %>%
        as.data.frame() %>%
        rownames_to_column("metric") %>%
        tibble() %>%
        mutate(pathogen = poi, .before = metric) %>%
        mutate(gene = goi, .before = metric) %>%
        mutate(Item = metric, .before = metric) %>%
        mutate(singular = ifelse(isSingular(model_abundances[[goi]][[poi]]), TRUE, FALSE)) %>%
        rbind(coeffs$abundances)
    }
  }
}


#####
# PATHOGEN PRESENCE VS TPM
pathogens_of_interest <- c("ABPV", "N. ceranae", "CBPV")

test_tibble_pathogen_presence <- list()
model_pathogen_presence <- list()
coeffs$pathogen_presence <- tibble()
for (goi in genes_of_interest) {
  for (poi in pathogens_of_interest) {
    # goi = "phosphoadenosine phosphosulfate reductase"
    # poi = "ABPV"
    test_tibble_pathogen_presence[[goi]][[poi]] <- gene_tpm %>%
      left_join(., metadata[c("Sample_ID", "Bee_pool", "Country", "Hive_ID", "Season", "Gut_part")], by = "Sample_ID") %>%
      left_join(., pathogen_ct, by = "Bee_pool", relationship = "many-to-many") %>%
      filter(gene == goi,
             pathogen == poi) %>%
      mutate(log_tpm = log10(tpm)) %>%
      filter(!is.infinite(log_tpm)) %>%
      mutate(pathogen_presence = ifelse(Ct < 41, 1, 0))
    
    model_pathogen_presence[[goi]][[poi]] <- glmer(
      pathogen_presence ~ log_tpm + Season + Gut_part + ( 1 | Hive_ID ),
      data = test_tibble_pathogen_presence[[goi]][[poi]],
      family = binomial)
    
    has_convergence_issues <- FALSE
    messages <- model_pathogen_presence[[goi]][[poi]]@optinfo$conv$lme4$messages
    if (!is.null(messages) & any(str_detect(messages, "failed to converge"))) {
      has_convergence_issues <- TRUE
    }
    if (has_convergence_issues) {
      print(paste0(goi," - ", poi, " - Model didn't converge. Removed from the list."))
      model_pathogen_presence[[goi]][[poi]] <- NULL
    } else {
      coeffs$pathogen_presence <- summary(model_pathogen_presence[[goi]][[poi]])$coefficients %>%
        as.data.frame() %>%
        rownames_to_column("metric") %>%
        tibble() %>%
        mutate(pathogen = poi, .before = metric) %>%
        mutate(gene = goi, .before = metric) %>%
        mutate(Item = metric, .before = metric) %>%
        mutate(singular = ifelse(isSingular(model_pathogen_presence[[goi]][[poi]]), TRUE, FALSE)) %>%
        rbind(coeffs$pathogen_presence)
    }
  }
}

#####
# PATHOGEN PRESENCE VS GENE PRESENCE
pathogens_of_interest <- c("ABPV", "N. ceranae", "CBPV")

test_tibble_both_presence <- list()
model_both_presence <- list()
coeffs$both_presence <- tibble()
for (goi in genes_of_interest) {
  for (poi in pathogens_of_interest) {
    # goi = "phosphoadenosine phosphosulfate reductase"
    # poi = "ABPV"
    test_tibble_both_presence[[goi]][[poi]] <- gene_tpm %>%
      left_join(., metadata[c("Sample_ID", "Bee_pool", "Country", "Hive_ID", "Season", "Gut_part")], by = "Sample_ID") %>%
      left_join(., pathogen_ct, by = "Bee_pool", relationship = "many-to-many") %>%
      filter(gene == goi,
             pathogen == poi) %>%
      mutate(gene_presence = ifelse(tpm > 0 , 1, 0)) %>%
      mutate(pathogen_presence = ifelse(Ct < 41, 1, 0))
    
    model_both_presence[[goi]][[poi]] <- glmer(
      pathogen_presence ~ gene_presence + Season + Gut_part + ( 1 | Hive_ID ),
      data = test_tibble_both_presence[[goi]][[poi]],
      family = binomial)
    
    has_convergence_issues <- FALSE
    messages <- model_both_presence[[goi]][[poi]]@optinfo$conv$lme4$messages
    if (!is.null(messages) & any(str_detect(messages, "failed to converge"))) {
      has_convergence_issues <- TRUE
    }
    if (has_convergence_issues) {
      print(paste0(goi," - ", poi, " - Model didn't converge. Removed from the list."))
      model_both_presence[[goi]][[poi]] <- NULL
    } else {
      coeffs$both_presence <- summary(model_both_presence[[goi]][[poi]])$coefficients %>%
        as.data.frame() %>%
        rownames_to_column("metric") %>%
        tibble() %>%
        mutate(pathogen = poi, .before = metric) %>%
        mutate(gene = goi, .before = metric) %>%
        mutate(Item = metric, .before = metric) %>%
        mutate(singular = ifelse(isSingular(model_both_presence[[goi]][[poi]]), TRUE, FALSE)) %>%
        rbind(coeffs$both_presence)
    }
  }
}

#####
# EXTRACT SLOPES
slopes <- list()
for (level in names(coeffs)) {
  slopes[[level]] <- coeffs[[level]] %>%
    filter(metric %in% c("Ct", "gene_presence", "log_tpm")) %>%
    rename(raw_p_value = any_of(c("Pr(>|z|)", "Pr(>|t|)"))) %>%
    mutate(raw_p_significant = case_when(raw_p_value <= 0.001 ~ "***",
                                         raw_p_value <= 0.01 ~ "**",
                                         raw_p_value <= 0.05 ~ "*",
                                         raw_p_value <= 0.075 ~ ".",
                                         .default = "n.s."
    )) %>%
    mutate(test_name = paste0(gene, "; ", pathogen, "; ", level), .before = pathogen) %>%
    mutate(p_adjusted = p.adjust(raw_p_value, method = "BH")) %>%
    mutate(p_adjust_significant = case_when(p_adjusted <= 0.001 ~ "***",
                                            p_adjusted <= 0.01 ~ "**",
                                            p_adjusted <= 0.05 ~ "*",
                                            p_adjusted <= 0.075 ~ ".",
                                            .default = "n.s."))
}


all_slopes <- bind_rows(slopes)

