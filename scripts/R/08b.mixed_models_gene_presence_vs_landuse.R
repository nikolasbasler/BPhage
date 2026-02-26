library(lme4)
library(lmerTest)
library(ggrepel)
library(forcats)
library(tidyverse)
library(patchwork)
library(DHARMa)

source("scripts/R/helpers/mixed_helpers.R")

genes_of_interest <- c("PAPS reductase")

contigs_with.PAPS.reductase <- read_lines("output/R/gene_content/amg_curation/contigs_with.PAPS reductase.txt")

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

metadata <- readRDS("data/metadata.RDS") %>%
  mutate(Hive_ID = as.character(Hive_ID))
classification <- readRDS("data/classification.RDS")

phage_tpm <- read.csv("output/R/relative_abundance/phage_tpm.csv") %>%
  tibble()

phold_predictions_with_extensions <- read.delim("output/R/gene_content/phold_predictions_with_extensions_bphage_renamed_genes.tsv")

gene_tpm <- phage_tpm %>%
  pivot_longer(-contig, names_to = "Sample_ID", values_to = "tpm") %>%
  mutate(tpm = ifelse(contig %in% contigs_with.PAPS.reductase, tpm, 0)) %>%
  group_by(Sample_ID) %>%
  mutate(tpm = sum(tpm),
         gene = "PAPS reductase") %>%
  ungroup() %>%
  select(Sample_ID, gene, tpm) %>%
  distinct()

coeffs_logit <- list()

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

model_logit <- list()
###### 
# CROPLAND

for (goi in genes_of_interest) {
  
  model_logit$Cropland_in_2km_radius[[goi]] <- glmer(
    presence ~ Cropland_in_2km_radius + Gut_part + Season + (1 | Hive_ID ),
    data = test_tibble_logit[[goi]],
    family = binomial)
  
  has_convergence_issues <- FALSE
  messages <- model_logit$Cropland_in_2km_radius[[goi]]@optinfo$conv$lme4$messages
  if (!is.null(messages) & any(str_detect(messages, "failed to converge"))) {
    has_convergence_issues <- TRUE
  }
  
  if (has_convergence_issues) {
    print(paste0(goi, " - Model didn't converge. Removed from the list."))
    model_logit$Cropland_in_2km_radius[[goi]] <- NULL
  } else {
    coeffs_logit$cropland <- summary(model_logit$Cropland_in_2km_radius[[goi]])$coefficients %>%
      as.data.frame() %>%
      rownames_to_column("metric") %>%
      tibble() %>%
      mutate(gene = goi, .before = metric) %>%
      mutate(Item = metric, .before = metric) %>%
      mutate(singular = ifelse(isSingular(model_logit$Cropland_in_2km_radius[[goi]]), TRUE, FALSE)) %>%
      rbind(coeffs_logit$cropland)
  }
}

######

##### 
# PESTICIDES:

coeffs_logit$total_pest <- tibble()
for (goi in genes_of_interest) {
  temp_test_tibble <- test_tibble_logit[[goi]] %>%
    rename(est_use_in_2k_radius = `Pesticides (total)`)
  
  model_logit$`Pesticides (total)`[[goi]] <- glmer(presence ~ est_use_in_2k_radius + Gut_part + Season +
                                                                   (1 | Hive_ID ), data = temp_test_tibble,
                                                                 family = binomial)
  
  has_convergence_issues <- FALSE
  messages <- model_logit$`Pesticides (total)`[[goi]]@optinfo$conv$lme4$messages
  if (!is.null(messages) & any(str_detect(messages, "failed to converge"))) {
    has_convergence_issues <- TRUE
  }
  
  if (has_convergence_issues) {
    print(paste0(goi, " - Model didn't converge. Removed from the list."))
    model_logit$`Pesticides (total)`[[goi]] <- NULL
  } else {
    coeffs_logit$total_pest <- summary(model_logit$`Pesticides (total)`[[goi]])$coefficients %>%
      as.data.frame() %>%
      rownames_to_column("metric") %>%
      tibble() %>%
      mutate(gene = goi, 
             Item = "Pesticides (total)",
             .before = metric) %>%
      mutate(singular = ifelse(isSingular(model_logit$`Pesticides (total)`[[goi]]), TRUE, FALSE)) %>%
      rbind(coeffs_logit$total_pest)
  }
}

#####
# PEST GROUPS
coeffs_logit$pest_groups <- tibble()
for (goi in genes_of_interest) {
  for (item in c("Insecticides", "Herbicides", "Fungicides and Bactericides", "Plant Growth Regulators")) {
    temp_test_tibble <- test_tibble_logit[[goi]] %>% 
      rename(est_use_in_2k_radius = all_of(item))
    
    model_logit[[item]][[goi]] <- glmer(presence ~ est_use_in_2k_radius + Gut_part + Season +
                                                      (1 | Hive_ID ), data = temp_test_tibble,
                                                    family = binomial)
    
    has_convergence_issues <- FALSE
    messages <- model_logit[[item]][[goi]]@optinfo$conv$lme4$messages
    if (!is.null(messages) & any(str_detect(messages, "failed to converge"))) {
      has_convergence_issues <- TRUE
    }
    
    if (has_convergence_issues) {
      print(paste0(goi, "; ", item, " - Model didn't converge. Removed from the list."))
      model_logit[[item]][[goi]] <- NULL
    } else {
      coeffs_logit$pest_groups <- summary(model_logit[[item]][[goi]])$coefficients %>%
        as.data.frame() %>%
        rownames_to_column("metric") %>%
        tibble() %>%
        mutate(gene = goi, 
               Item = item,
               .before = metric) %>%
        mutate(singular = ifelse(isSingular(model_logit[[item]][[goi]]), TRUE, FALSE)) %>%
        rbind(coeffs_logit$pest_groups)
    }
  }
}

#####
# SPECIFIC PESTS

spec_pests <- tibble(Item = colnames(cropland_and_FAO)) %>%
  filter(str_detect(Item, "Herbicides ") | str_detect(Item, "Fung & Bact ") | str_detect(Item, "Insecticides ")) %>%
  unlist(use.names = FALSE)

coeffs_logit$specific_pests <- tibble()
for (goi in genes_of_interest) {
  for (item in spec_pests) {
    temp_test_tibble <- test_tibble_logit[[goi]] %>%
      rename(est_use_in_2k_radius = all_of(item))
    
    model_logit[[item]][[goi]] <- glmer(presence ~ est_use_in_2k_radius + Gut_part + Season +
                                                         (1 | Hive_ID ), data = temp_test_tibble,
                                                       family = binomial)
    
    has_convergence_issues <- FALSE
    messages <- model_logit[[item]][[goi]]@optinfo$conv$lme4$messages
    if (!is.null(messages) & any(str_detect(messages, "failed to converge"))) {
      has_convergence_issues <- TRUE
    }
    if (has_convergence_issues) {
      print(paste0(goi, "; ", item, " - Model didn't converge. Removed from the list."))
      model_logit[[item]][[goi]] <- NULL
    } else {
      coeffs_logit$specific_pests <- summary(model_logit[[item]][[goi]])$coefficients %>%
        as.data.frame() %>%
        rownames_to_column("metric") %>%
        tibble() %>%
        mutate(gene = goi,
               Item = item,
               .before = metric) %>%
        mutate(singular = ifelse(isSingular(model_logit[[item]][[goi]]), TRUE, FALSE)) %>%
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
             se_intercept = `Std. Error`) %>%
      select(test_name, intercept, se_intercept) %>%
      full_join(temp_slope_tibble, ., by = "test_name") %>%
      relocate(c(intercept, se_intercept), .after = `Std. Error`)
    
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

color_list <- list(dark = list(`PAPS reductase` = "#ef8f01"),
                   bright = list(`PAPS reductase` = "#ef8f01"))

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


# All results
all_tests_forest_plot <- all_slopes %>%
  mutate(axis_labels = fct_rev(fct_inorder(test_name))) %>%
  mutate(estimate = Estimate,
         error = `Std. Error`) %>%
  forest_plot(plot_title = "all tests")

#####
# DIAGNOSTICS

model_diagnostics <- list()
for (item in names(model_logit)) {
  for (goi in names(model_logit[[item]])) {
    slope_of_interest <- all_slopes %>% 
      filter(Item == item,
             gene == goi) %>%
      reframe(test_name = test_name,
              is_significant = p_adjusted < 0.05)
    if (slope_of_interest$is_significant) {
      model_diagnostics[[item]][[goi]] <- diagnostics_logistic_model(model_logit[[item]][[goi]], name = slope_of_interest$test_name)
    }
  }
}

#####
# SAVE FILES

system("mkdir -p output/R/genes_pathogens_and_landuse/gene_presence_vs_landuse/single_panels")
system("mkdir -p output/R/genes_pathogens_and_landuse/gene_presence_vs_landuse/model_diagnostics")

write_delim(cropland_and_FAO, "output/R/genes_pathogens_and_landuse/cropland_and_FAO.tsv",
            delim = "\t")

write_delim(bind_rows(coeffs_logit), "output/R/genes_pathogens_and_landuse/gene_presence_vs_landuse/gene_presence_vs_landuse.all_coeffs.tsv",
            delim = "\t")
write_delim(all_slopes, "output/R/genes_pathogens_and_landuse/gene_presence_vs_landuse/gene_presence_vs_landuse.all_slopes.tsv",
            delim = "\t")
ggsave("output/R/genes_pathogens_and_landuse/gene_presence_vs_landuse/gene_presence_vs_landuse.all_tests.pdf",
       all_tests_forest_plot, width = 12, height = 40, limitsize = FALSE)

for (goi in names(genes_logit_plots)) {
  for (item in names(genes_logit_plots[[goi]])) {
    ggsave(paste0("output/R/genes_pathogens_and_landuse/gene_presence_vs_landuse/single_panels/", goi, ".", item, ".pdf"),
           genes_logit_plots[[goi]][[item]], height = 3.5, width = 3.5)
    saveRDS(genes_logit_plots[[goi]][[item]],
            paste0("output/R/genes_pathogens_and_landuse/gene_presence_vs_landuse/single_panels/RDS.", goi, ".", item, ".rds"))
  }
}
ggsave("output/R/genes_pathogens_and_landuse/gene_presence_vs_landuse/single_panels/common_legend.pdf",
       common_legend, height = 1, width = 6)

for (item in names(model_diagnostics)) {
  for (goi in names(model_diagnostics[[item]])) {
    ggsave(paste0("output/R/genes_pathogens_and_landuse/gene_presence_vs_landuse/model_diagnostics/model_diagnostics.", item, ".", goi, ".pdf"),
           model_diagnostics[[item]][[goi]], width = 6, height = 6)
  }
}

