library(lme4)
library(lmerTest)
library(ggrepel)
library(forcats)
library(tidyverse)
library(patchwork)
library(tidytext)

source("scripts/R/helpers/mixed_helpers.R")

cropland_fraction <- read.csv("data/land_cover_results.csv") %>% 
  tibble() %>%
  mutate(cropland_fraction = cropland_fraction / 100) %>%
  rename(cropland_fraction_2k_radius = cropland_fraction) %>%
  arrange(cropland_fraction_2k_radius)
# make FAOSTAT_added_data:
source("scripts/R/helpers/FAOstat_table.R")

nosema_relabund <- read.delim("output/nosema_mapped_counts_all.tsv") %>% 
  tibble() %>%
  mutate(tpm = mapped_to_Nosema / Hostout_R1_plus_R2) %>%
  select(Sample_ID, tpm)

cropland_and_FAO <- FAOSTAT_added_data %>%
  select(Country, Item, est_use_in_2k_radius) %>%
  pivot_wider(id_cols = Country, values_from = est_use_in_2k_radius, names_from = Item) %>%
  # left_join(cropland_fraction[c("Country", "cropland_fraction_2k_radius")], ., by = "Country") %>%
  left_join(FAOSTAT_added_data[c("Country", "cropland_fraction_2k_radius", "Cropland_in_2km_radius")],. , by = "Country") %>%
  distinct() %>%
  select_if(~ !any(is.na(.)))

metadata <- readRDS("output/R/R_variables/metadata.RDS") %>%
  mutate(Hive_ID = as.character(Hive_ID))
classification <- readRDS("output/R/R_variables/classification.RDS")

phage_tpm <- read.csv("output/R/relative_abundance/phage_tpm.csv") %>%
  tibble()

phold_predictions_with_extensions <- read.delim("output/R/gene_content/phold_predictions_with_extensions_bphage_renamed_genes.tsv")

kegg_mapping <- read.delim("data/kegg_mapping.tsv", colClasses = "character") %>%
  tibble()

kegg_and_phold <- kegg_mapping %>%
  left_join(., phold_predictions_with_extensions[c("cds_id", "phrog", "function.", "product")], by = "cds_id")

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

genes_of_interest <- c("nosema_relabund")

coeffs_tpm_simple <- list()

model_tpm <- list()

#####
# CROPLAND
test_tibble_log_tpm <- list()
coeffs_tpm_simple$cropland <- tibble()
for (goi in genes_of_interest) {
  
  test_tibble_log_tpm[[goi]] <- nosema_relabund %>%
    # filter(gene == goi) %>%
    left_join(., metadata[c("Sample_ID", "Country", "Hive_ID", "Season", "Gut_part")], by = "Sample_ID") %>%
    distinct() %>%
    left_join(., cropland_and_FAO, by = "Country") %>%
    mutate(Gut_part = factor(Gut_part, levels = c("rec", "ile", "mid"))) %>%
    mutate(log_tpm = log10(tpm)) %>%
    filter(!is.infinite(log_tpm))
  
  model_tpm$Cropland_in_2km_radius[[goi]] <- lmer(log_tpm ~ Cropland_in_2km_radius + Gut_part + Season +
                                                    (1 | Hive_ID ), data = test_tibble_log_tpm[[goi]])
  
  has_convergence_issues <- FALSE
  messages <- model_tpm$Cropland_in_2km_radius[[goi]]@optinfo$conv$lme4$messages
  if (is.null(messages)) {
    has_convergence_issues <- FALSE
  } else {
    has_convergence_issues <- length(messages[!grepl("singular", messages, ignore.case = TRUE)]) != 0
  }
  if (has_convergence_issues) {
    print(paste0(goi, " - Simple model didn't converge. Removed from the list."))
    model_tpm$Cropland_in_2km_radius[[goi]] <- NULL
    
  } else {
    coeffs_tpm_simple$cropland <- summary(model_tpm$Cropland_in_2km_radius[[goi]])$coefficients %>%
      as.data.frame() %>%
      rownames_to_column("metric") %>%
      tibble() %>%
      mutate(gene = goi, .before = metric) %>%
      mutate(Item = metric, .before = metric) %>%
      mutate(singular = ifelse(isSingular(model_tpm$Cropland_in_2km_radius[[goi]]), TRUE, FALSE)) %>%
      rbind(coeffs_tpm_simple$cropland)
  }
}

##### 
# PESTICIDES:

coeffs_tpm_simple$total_pest <- tibble()
for (goi in genes_of_interest) {
  temp_test_tibble <- test_tibble_log_tpm[[goi]] %>% 
    rename(est_use_in_2k_radius = "Pesticides (total)")
  model_tpm$`Pesticides (total)`[[goi]] <- lmer(log_tpm ~ est_use_in_2k_radius + Gut_part + Season +
                                                  (1 | Hive_ID ), data = temp_test_tibble)
  
  has_convergence_issues <- FALSE
  messages <- model_tpm$`Pesticides (total)`[[goi]]@optinfo$conv$lme4$messages
  if (is.null(messages)) {
    has_convergence_issues <- FALSE
  } else {
    has_convergence_issues <- length(messages[!grepl("singular", messages, ignore.case = TRUE)]) != 0
  }
  if (has_convergence_issues) {
    print(paste0(goi, " - Simple model didn't converge. Removed from the list."))
    model_tpm$`Pesticides (total)`[[goi]] <- NULL
    
  } else {
    coeffs_tpm_simple$total_pest <- summary(model_tpm$`Pesticides (total)`[[goi]])$coefficients %>%
      as.data.frame() %>%
      rownames_to_column("metric") %>%
      tibble() %>%
      mutate(gene = goi, 
             Item = "Pesticides (total)",
             .before = metric) %>%
      mutate(singular = ifelse(isSingular(model_tpm$`Pesticides (total)`[[goi]]), TRUE, FALSE)) %>%
      rbind(coeffs_tpm_simple$total_pest)
  }
}


#####
# PEST GROUPS
coeffs_tpm_simple$pest_groups <- tibble()
for (goi in genes_of_interest) {
  for (item in c("Insecticides", "Herbicides", "Fungicides and Bactericides", "Plant Growth Regulators")) {
    temp_test_tibble <- test_tibble_log_tpm[[goi]] %>% 
      rename(est_use_in_2k_radius = all_of(item))
    model_tpm[[item]][[goi]] <- lmer(log_tpm ~ est_use_in_2k_radius + Gut_part + Season +
                                       (1 | Hive_ID ), data = temp_test_tibble)
    
    has_convergence_issues <- FALSE
    messages <- model_tpm[[item]][[goi]]@optinfo$conv$lme4$messages
    if (is.null(messages)) {
      has_convergence_issues <- FALSE
    } else {
      has_convergence_issues <- length(messages[!grepl("singular", messages, ignore.case = TRUE)]) != 0
    }
    if (has_convergence_issues) {
      print(paste0(goi, " - Simple model didn't converge. Removed from the list."))
      model_tpm[[item]][[goi]] <- NULL
      
    } else {
      coeffs_tpm_simple$pest_groups <- summary(model_tpm[[item]][[goi]])$coefficients %>%
        as.data.frame() %>%
        rownames_to_column("metric") %>%
        tibble() %>%
        mutate(gene = goi, 
               Item = item,
               .before = metric) %>%
        mutate(singular = ifelse(isSingular(model_tpm[[item]][[goi]]), TRUE, FALSE)) %>%
        rbind(coeffs_tpm_simple$pest_groups)
    }
  }
}

#####
# SPECIFIC PESTS

spec_pests <- tibble(Item = colnames(cropland_and_FAO)) %>%
  filter(str_detect(Item, "Herbicides ") | str_detect(Item, "Fung & Bact ") | str_detect(Item, "Insecticides ")) %>%
  unlist(use.names = FALSE)

model_tpm_simple_specific_pests <-list()
coeffs_tpm_simple$specific_pests <- tibble()
for (goi in genes_of_interest) {
  for (item in spec_pests) {
    temp_test_tibble <- test_tibble_log_tpm[[goi]] %>% 
      rename(est_use_in_2k_radius = all_of(item))
    model_tpm[[item]][[goi]] <- lmer(log_tpm ~ est_use_in_2k_radius + Gut_part + Season +
                                       (1 | Hive_ID ), data = temp_test_tibble)
    
    has_convergence_issues <- FALSE
    messages <- model_tpm[[item]][[goi]]@optinfo$conv$lme4$messages
    if (is.null(messages)) {
      has_convergence_issues <- FALSE
    } else {
      has_convergence_issues <- length(messages[!grepl("singular", messages, ignore.case = TRUE)]) != 0
    }
    if (has_convergence_issues) {
      print(paste0(goi, " - Simple model didn't converge. Removed from the list."))
      model_tpm[[item]][[goi]] <- NULL
      
    } else {
      coeffs_tpm_simple$specific_pests <- summary(model_tpm[[item]][[goi]])$coefficients %>%
        as.data.frame() %>%
        rownames_to_column("metric") %>%
        tibble() %>%
        mutate(gene = goi, 
               Item = item,
               .before = metric) %>%
        mutate(singular = ifelse(isSingular(model_tpm[[item]][[goi]]), TRUE, FALSE)) %>%
        rbind(coeffs_tpm_simple$specific_pests)
    }
    
  }
}

#####
# EXTRACT SLOPES
slopes <- list()
for (level in names(coeffs_tpm_simple)) {
  
  temp_slope_tibble <- coeffs_tpm_simple[[level]] %>%
    filter(metric %in% c("Cropland_in_2km_radius", "est_use_in_2k_radius")) %>%
    rename(raw_p_value = `Pr(>|t|)`) %>%
    mutate(raw_p_significant = case_when(raw_p_value <= 0.001 ~ "***",
                                         raw_p_value <= 0.01 ~ "**",
                                         raw_p_value <= 0.05 ~ "*",
                                         raw_p_value <= 0.075 ~ ".",
                                         .default = "n.s."
    )) %>%
    mutate(test_name = paste0(gene, "; ", Item), .before = gene)
  
  slopes[[level]] <- coeffs_tpm_simple[[level]] %>%
    filter(metric == "(Intercept)") %>%
    mutate(test_name = paste0(gene, "; ", temp_slope_tibble$Item), .before = gene) %>%
    rename(intercept = Estimate,
           sd_intercept = `Std. Error`) %>%
    select(test_name, intercept, sd_intercept) %>%
    full_join(temp_slope_tibble, ., by = "test_name") %>%
    relocate(c(intercept, sd_intercept), .after = `Std. Error`)
  
}

##### 
# ADJUST P-VALUES

all_slopes <- bind_rows(slopes) %>%
  mutate(p_adjusted = p.adjust(raw_p_value, method = "BH")) %>%
  mutate(p_adjust_significant = case_when(p_adjusted <= 0.001 ~ "***",
                                          p_adjusted <= 0.01 ~ "**",
                                          p_adjusted <= 0.05 ~ "*",
                                          p_adjusted <= 0.075 ~ ".",
                                          .default = "n.s.")
  )


#####
# DIAGNOSTICS

model_diagnostics <- list()
for (item in names(model_tpm)) {
  for (goi in names(model_tpm[[item]])) {
    slope_of_interest <- all_slopes %>%
      filter(Item == item,
             gene == goi) %>%
      reframe(test_name = test_name,
              is_significant = p_adjusted < 0.05)
    if (slope_of_interest$is_significant) {
      model_diagnostics[[item]][[goi]] <- diagnostics_linear_model(model_tpm[[item]][[goi]], name = slope_of_interest$test_name)
    }
  }
}

# All results
all_tests_forest_plot <- all_slopes %>%
  mutate(axis_labels = fct_rev(fct_inorder(test_name))) %>%
  mutate(estimate = Estimate,
         error = `Std. Error`) %>%
  forest_plot(plot_title = "all tests")



###
# Only one spec pest is significant after BH correction (negative effect), and the diagnostics don't look too healthy:
all_tests_forest_plot
all_slopes %>%
  filter(Item == "Herbicides – Dinitroanilines")
model_diagnostics
