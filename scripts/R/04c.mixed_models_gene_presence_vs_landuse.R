library(lme4)
library(lmerTest)
library(ggrepel)
library(forcats)
library(tidyverse)
library(patchwork)
library(DHARMa)

source("scripts/R/helpers/mixed_helpers.R")

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

genes_of_interest <- c("Chitinase",
                       "Glucosyltransferase",
                       "Levanase",
                       "PAPS reductase", 
                       "PnuC")

coeffs_logit <- list()
countries_with_presence <- list()
samples_with_presence <- list()

###### 
# TEST TIBBLE
test_tibble_logit <- list()
for (goi in genes_of_interest) {
  
  test_tibble_logit[[goi]] <- gene_tpm %>%
    filter(gene == goi) %>%
    left_join(., metadata[c("Sample_ID", "Country", "Hive_ID", "Season", "Gut_part")], by = "Sample_ID") %>%
    distinct() %>%
    left_join(., cropland_and_FAO, by = "Country") %>%
    mutate(Gut_part = factor(Gut_part, levels = c("rec", "ile", "mid"))) %>%
    mutate(presence = ifelse(tpm > 0, 1, 0), .before = tpm)
}

###### 
# CROPLAND

model_logit_cropland <- list()
samples_with_presence$cropland <- tibble()
for (goi in genes_of_interest) {
  
  model_logit_cropland[[goi]] <- glmer(
    presence ~ Cropland_in_2km_radius + Gut_part + Season + (1 | Hive_ID ),
    data = test_tibble_logit[[goi]],
    family = binomial)
  
  has_convergence_issues <- FALSE
  messages <- model_logit_cropland[[goi]]@optinfo$conv$lme4$messages
  if (!is.null(messages) & any(str_detect(messages, "failed to converge"))) {
    has_convergence_issues <- TRUE
  }
  
  if (has_convergence_issues) {
    print(paste0(goi, " - Model didn't converge. Removed from the list."))
    model_logit_cropland[[goi]] <- NULL
  } else {
    coeffs_logit$cropland <- summary(model_logit_cropland[[goi]])$coefficients %>%
      as.data.frame() %>%
      rownames_to_column("metric") %>%
      tibble() %>%
      mutate(gene = goi, .before = metric) %>%
      mutate(Item = metric, .before = metric) %>%
      mutate(singular = ifelse(isSingular(model_logit_cropland[[goi]]), TRUE, FALSE)) %>%
      rbind(coeffs_logit$cropland)
  }
}

######

##### 
# PESTICIDES:

model_logit_total_pest <- list()
coeffs_logit$total_pest <- tibble()
for (goi in genes_of_interest) {
  temp_test_tibble <- test_tibble_logit[[goi]] %>%
    rename(est_use_in_2k_radius = `Pesticides (total)`)
  
  model_logit_total_pest[[goi]][["Pesticides (total)"]] <- glmer(presence ~ est_use_in_2k_radius + Gut_part + Season +
                                                                   (1 | Hive_ID ), data = temp_test_tibble,
                                                                 family = binomial)
  
  has_convergence_issues <- FALSE
  messages <- model_logit_total_pest[[goi]][["Pesticides (total)"]]@optinfo$conv$lme4$messages
  if (!is.null(messages) & any(str_detect(messages, "failed to converge"))) {
    has_convergence_issues <- TRUE
  }
  
  if (has_convergence_issues) {
    print(paste0(goi, " - Model didn't converge. Removed from the list."))
    model_logit_total_pest[[goi]] <- NULL
  } else {
    coeffs_logit$total_pest <- summary(model_logit_total_pest[[goi]][["Pesticides (total)"]])$coefficients %>%
      as.data.frame() %>%
      rownames_to_column("metric") %>%
      tibble() %>%
      mutate(gene = goi, 
             Item = "Pesticides (total)",
             .before = metric) %>%
      mutate(singular = ifelse(isSingular(model_logit_total_pest[[goi]][["Pesticides (total)"]]), TRUE, FALSE)) %>%
      rbind(coeffs_logit$total_pest)
  }
}


#####
# PEST GROUPS
model_logit_pest_groups <- list()
coeffs_logit$pest_groups <- tibble()
for (goi in genes_of_interest) {
  for (item in c("Insecticides", "Herbicides", "Fungicides and Bactericides", "Plant Growth Regulators")) {
    temp_test_tibble <- test_tibble_logit[[goi]] %>% 
      rename(est_use_in_2k_radius = all_of(item))
    
    model_logit_pest_groups[[goi]][[item]] <- glmer(presence ~ est_use_in_2k_radius + Gut_part + Season +
                                                      (1 | Hive_ID ), data = temp_test_tibble,
                                                    family = binomial)
    
    has_convergence_issues <- FALSE
    messages <- model_logit_pest_groups[[goi]][[item]]@optinfo$conv$lme4$messages
    if (!is.null(messages) & any(str_detect(messages, "failed to converge"))) {
      has_convergence_issues <- TRUE
    }
    
    if (has_convergence_issues) {
      print(paste0(goi, "; ", item, " - Model didn't converge. Removed from the list."))
      model_logit_pest_groups[[goi]][[item]] <- NULL
    } else {
      coeffs_logit$pest_groups <- summary(model_logit_pest_groups[[goi]][[item]])$coefficients %>%
        as.data.frame() %>%
        rownames_to_column("metric") %>%
        tibble() %>%
        mutate(gene = goi, 
               Item = item,
               .before = metric) %>%
        mutate(singular = ifelse(isSingular(model_logit_pest_groups[[goi]][[item]]), TRUE, FALSE)) %>%
        rbind(coeffs_logit$pest_groups)
    }
  }
}

#####
# SPECIFIC PESTS

spec_pests <- tibble(Item = colnames(cropland_and_FAO)) %>%
  filter(str_detect(Item, "Herbicides ") | str_detect(Item, "Fung & Bact ") | str_detect(Item, "Insecticides ")) %>%
  unlist(use.names = FALSE)

model_logit_specific_pests <- list()
coeffs_logit$specific_pests <- tibble()
for (goi in genes_of_interest) {
  for (item in spec_pests) {
    temp_test_tibble <- test_tibble_logit[[goi]] %>%
      rename(est_use_in_2k_radius = all_of(item))
    
    model_logit_specific_pests[[goi]][[item]] <- glmer(presence ~ est_use_in_2k_radius + Gut_part + Season +
                                                         (1 | Hive_ID ), data = temp_test_tibble,
                                                       family = binomial)
    
    has_convergence_issues <- FALSE
    messages <- model_logit_specific_pests[[goi]][[item]]@optinfo$conv$lme4$messages
    if (!is.null(messages) & any(str_detect(messages, "failed to converge"))) {
      has_convergence_issues <- TRUE
    }
    if (has_convergence_issues) {
      print(paste0(goi, "; ", item, " - Model didn't converge. Removed from the list."))
      model_logit_pest_groups[[goi]][[item]] <- NULL
    } else {
      coeffs_logit$specific_pests <- summary(model_logit_specific_pests[[goi]][[item]])$coefficients %>%
        as.data.frame() %>%
        rownames_to_column("metric") %>%
        tibble() %>%
        mutate(gene = goi,
               Item = item,
               .before = metric) %>%
        mutate(singular = ifelse(isSingular(model_logit_specific_pests[[goi]][[item]]), TRUE, FALSE)) %>%
        rbind(coeffs_logit$specific_pests)
    }
  }
}

#####
# EXTRACT SLOPES
slopes <- list()
for (level in names(coeffs_logit)) {
  if (nrow(coeffs_logit[[level]]) > 0) {
    
    temp_slope_tibble <- coeffs_logit[[level]] %>%
      filter(metric %in% c("Cropland_in_2km_radius", "est_use_in_2k_radius")) %>%
      rename(raw_p_value = `Pr(>|z|)`) %>%
      mutate(raw_p_significant = case_when(raw_p_value <= 0.001 ~ "***",
                                           raw_p_value <= 0.01 ~ "**",
                                           raw_p_value <= 0.05 ~ "*",
                                           raw_p_value <= 0.075 ~ ".",
                                           .default = "n.s."
      )) %>%
      mutate(test_name = paste0(gene, "; ", Item), .before = gene)
    
    slopes[[level]] <- coeffs_logit[[level]] %>%
      filter(metric == "(Intercept)") %>%
      mutate(test_name = paste0(gene, "; ", temp_slope_tibble$Item), .before = gene) %>%
      rename(intercept = Estimate,
             sd_intercept = `Std. Error`) %>%
      select(test_name, intercept, sd_intercept) %>%
      full_join(temp_slope_tibble, ., by = "test_name") %>%
      relocate(c(intercept, sd_intercept), .after = `Std. Error`)
    
  }
}

all_slopes <- bind_rows(slopes) %>%
  mutate(p_adjusted = p.adjust(raw_p_value, method = "BH")) %>%
  mutate(p_adjust_significant = case_when(p_adjusted <= 0.001 ~ "***",
                                          p_adjusted <= 0.01 ~ "**",
                                          p_adjusted <= 0.05 ~ "*",
                                          p_adjusted <= 0.075 ~ ".",
                                          .default = "n.s."))


#####
# MAKE PLOTS

# Significant results
lowest_highest <- cropland_and_FAO %>%
  pivot_longer(-Country, names_to = "Item") %>%
  group_by(Item) %>%
  summarise(lowest = min(value),
            highest = max(value))

sig_tests <- all_slopes %>%
  filter(p_adjusted < 0.05) %>%
  left_join(., lowest_highest, by = "Item") %>%
  mutate(effect = logistic_effect_fun(s = Estimate, h = highest, l = lowest),
         effect = ifelse(effect < 1, 1 / effect, effect)) %>%
  group_by(gene) %>%
  mutate(y_stretching_factor = max(effect) / effect) %>%
  mutate(which_y_end_to_stretch = if_else(
    Estimate[y_stretching_factor == 1] > 0,
    "upper_end",
    "lower_end")) %>% 
  ungroup() %>%
  mutate(y_stretching_factor = case_when(
    which_y_end_to_stretch == "upper_end" & Estimate < 0 ~ 1/y_stretching_factor,
    which_y_end_to_stretch == "lower_end" & Estimate > 0 ~ 1/y_stretching_factor,
    .default = y_stretching_factor
  )) %>%
  select(-effect) 

color_list <- list(dark = list(Chitinase = "#8B4513", Glucosyltransferase = "#FFC300", `PAPS reductase` = "#ef8f01"),
                   bright = list(Chitinase = "#8B4513", Glucosyltransferase = "#FFC300", `PAPS reductase` = "#ef8f01"))

genes_logit_plots <- list()
for (t_name in unique(sig_tests$test_name)) {

  logit_test <- sig_tests %>%
    filter(test_name == t_name)
  
  tested_gene <- logit_test$gene
  tested_item <- logit_test$Item
  
  genes_logit_plots[[tested_gene]][[tested_item]] <- 
    mixed_model_plot(filt_test_tibble = logit_test,
                     transform_fun = logistic_fun,
                     effect_fun = logistic_effect_fun,
                     dark_col = color_list$dark[[tested_gene]],
                     bright_col = color_list$bright[[tested_gene]],
                     y_axis_label = "Probability")
}


common_legend <- legend_factory(title = "Gene", 
                      items = names(color_list$dark),
                      colors = unlist(color_list$dark),
                      position = "bottom")

wrap_of_wraps <- wrap_plots(
  wrap_plots(genes_logit_plots$Chitinase[1:3], nrow = 1, axes = "collect"),
    wrap_plots(genes_logit_plots$Chitinase[4:6], nrow = 1, axes = "collect"),
    wrap_plots(genes_logit_plots$Chitinase[7:9], nrow = 1, axes = "collect"),
    wrap_plots(list(genes_logit_plots$Chitinase$`Insecticides – Pyrethroids`,
                    genes_logit_plots$Chitinase$`Insecticides – Carbamates`, 
                    genes_logit_plots$Glucosyltransferase$`Insecticides – Organo-phosphates`), nrow = 1, axes = "collect"),
    wrap_plots(genes_logit_plots$`PAPS reductase`, nrow = 1, axes = "collect"),
    common_legend,
  nrow = 6, heights = c(rep(4, 5), 1)
)

# All results
all_tests_forest_plot <- all_slopes %>%
  mutate(axis_labels = fct_rev(fct_inorder(test_name))) %>%
  mutate(estimate = Estimate,
         error = `Std. Error`) %>%
  forest_plot(plot_title = "all tests")

#####
# DIAGNISTICS


model_diagnostics <- list()
for (goi in genes_of_interest) {
  model_diagnostics$cropland[[goi]] <- diagnostics_logistic_model(model_logit_cropland[[goi]], name = paste0(goi, "presence vs. Cropland_in_2km_radius"))
  model_diagnostics$total_pest[[goi]] <- diagnostics_logistic_model(model_logit_total_pest[[goi]], name = paste0(goi, "presence vs. Pesticides (total)"))
  for (computed_model in names(model_logit_pest_groups)) {
    model_logit_pest_groups
    
  }
  
  for (computed_model in names(model_logit_specific_pests)) {
    model_logit_specific_pests
    
  }
}


diagnostics_logistic_model

#####
# SAVE FILES

system("mkdir -p output/R/genes_pathogens_and_landuse/gene_presence_vs_landuse/single_panels")

write_delim(bind_rows(coeffs_logit), "output/R/genes_pathogens_and_landuse/gene_presence_vs_landuse/gene_presence_vs_landuse.all_coeffs.tsv",
            delim = "\t")
write_delim(all_slopes, "output/R/genes_pathogens_and_landuse/gene_presence_vs_landuse/gene_presence_vs_landuse.all_sloppes.tsv",
            delim = "\t")
ggsave("output/R/genes_pathogens_and_landuse/gene_presence_vs_landuse/gene_presence_vs_landuse.all_tests.pdf",
       all_tests_forest_plot, width = 12, height = 40, limitsize = FALSE)

ggsave("output/R/genes_pathogens_and_landuse/gene_presence_vs_landuse/gene_presence_vs_landuse.wrap.pdf",
       wrap_of_wraps, height = 13, width = 9.25)

for (goi in names(genes_logit_plots)) {
  for (item in names(genes_logit_plots[[goi]])) {
    ggsave(paste0("output/R/genes_pathogens_and_landuse/gene_presence_vs_landuse/single_panels/", goi, ".", item, ".pdf"),
           genes_logit_plots[[goi]][[item]], height = 3.5, width = 3.5)
  }
}
ggsave("output/R/genes_pathogens_and_landuse/gene_presence_vs_landuse/single_panels/common_legend.pdf",
       common_legend, height = 1, width = 6)



