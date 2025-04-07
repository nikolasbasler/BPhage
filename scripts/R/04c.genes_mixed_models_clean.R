library(lme4)
library(lmerTest)
library(DHARMa)
library(ggrepel)
library(forcats)
library(tidyverse)
library(patchwork)
library(gMCPLite)

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
  # left_join(cropland_fraction[c("Country", "cropland_fraction_2k_radius")], ., by = "Country") %>%
  left_join(FAOSTAT_added_data[c("Country", "cropland_fraction_2k_radius", "ha_cropland_in_2k_radius")],. , by = "Country") %>%
  distinct() %>%
  select_if(~ !any(is.na(.)))

metadata <- readRDS("output/R/R_variables/metadata.RDS") %>%
  mutate(Hive_ID = as.character(Hive_ID))
classification <- readRDS("output/R/R_variables/classification.RDS")

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

genes_of_interest <- unique(gene_tpm$gene)
coeffs <- list()

test_tibble_log_tpm <- list()
model_cropland <- list()
coeffs$cropland <- tibble()
for (goi in genes_of_interest) {
  
  test_tibble_log_tpm[[goi]] <- gene_tpm %>%
    filter(gene == goi) %>%
    left_join(., metadata[c("Sample_ID", "Country", "Hive_ID", "Season", "Gut_part")], by = "Sample_ID") %>%
    distinct() %>%
    left_join(., cropland_and_FAO, by = "Country") %>%
    mutate(Gut_part = factor(Gut_part, levels = c("rec", "ile", "mid"))) %>%
    mutate(log_tpm = log10(tpm)) %>%
    filter(!is.infinite(log_tpm)) %>%
    group_by(Hive_ID) %>%
    # filter(n() > 4) %>%
    ungroup()
  
  countries_left <- test_tibble_log_tpm[[goi]] %>%
    distinct(Country) %>%
    nrow()
  
  if (countries_left == 8 ) {
    model_cropland[[goi]] <- lmer(log_tpm ~ ha_cropland_in_2k_radius + Gut_part + Season +
                                   (1 | Hive_ID), data = test_tibble_log_tpm[[goi]])
    # summary(model_cropland[[goi]])
    
    # model_cropland[[goi]] <- lmer(log_tpm ~ ha_cropland_in_2k_radius + Gut_part + Season +
    #                                (1 | Hive_ID / Country),
    #                              data = test_tibble_log_tpm[[goi]])
    # summary(model_cropland[[goi]])
    coeffs$cropland <- summary(model_cropland[[goi]])$coefficients %>%
      as.data.frame() %>%
      rownames_to_column("metric") %>%
      tibble() %>%
      mutate(gene = goi, .before = metric) %>%
      mutate(Item = metric, .before = metric) %>%
      mutate(singular = ifelse(isSingular(model_cropland[[goi]]), TRUE, FALSE)) %>%
      rbind(coeffs$cropland)
  } else {
    print(paste0(goi, " - no model calculated. Only ", countries_left, " countries left."))
  }
}

genes_of_particular_interest <- names(model_cropland)

# summary(model_cropland$`phosphoadenosine phosphosulfate reductase`)
# summary(model_cropland$levanase)
# summary(model_cropland$glucosyltransferase)
# summary(model_cropland$chitinase)
# summary(model_cropland$`PnuC-like nicotinamide mononucleotide transport`)

##### 
# PESTICIDES:

total_pest_model <- list()
coeffs$total_pest <- tibble()
for (gopi in genes_of_particular_interest) {
  temp_test_tibble <- test_tibble_log_tpm[[gopi]] %>% 
    rename(est_use_in_2k_radius = "Pesticides (total)")
  total_pest_model[[gopi]][["Pesticides (total)"]] <- lmer(log_tpm ~ est_use_in_2k_radius + Gut_part + Season +
                                                     (1 | Hive_ID), data = temp_test_tibble)
  
  coeffs$total_pest <- summary(total_pest_model[[gopi]][["Pesticides (total)"]])$coefficients %>%
    as.data.frame() %>%
    rownames_to_column("metric") %>%
    tibble() %>%
    mutate(gene = gopi, 
           Item = "Pesticides (total)",
           .before = metric) %>%
    mutate(singular = ifelse(isSingular(total_pest_model[[gopi]][["Pesticides (total)"]]), TRUE, FALSE)) %>%
    rbind(coeffs$total_pest)
}

#####
# PEST GROUPS
pest_groups_model <- list()
coeffs$pest_groups <- tibble()
for (gopi in genes_of_particular_interest) {
  for (item in c("Insecticides", "Herbicides", "Fungicides and Bactericides", "Plant Growth Regulators")) {
    temp_test_tibble <- test_tibble_log_tpm[[gopi]] %>% 
      rename(est_use_in_2k_radius = all_of(item))
    pest_groups_model[[gopi]][[item]] <- lmer(log_tpm ~ est_use_in_2k_radius + Gut_part + Season +
                                                       (1 | Hive_ID), data = temp_test_tibble)
    
    coeffs$pest_groups <- summary(pest_groups_model[[gopi]][[item]])$coefficients %>%
      as.data.frame() %>%
      rownames_to_column("metric") %>%
      tibble() %>%
      mutate(gene = gopi, 
             Item = item,
             .before = metric) %>%
      mutate(singular = ifelse(isSingular(pest_groups_model[[gopi]][[item]]), TRUE, FALSE)) %>%
      rbind(coeffs$pest_groups)
  }
}


spec_pests <- tibble(Item = colnames(cropland_and_FAO)) %>%
  filter(str_detect(Item, "Herbicides ") | str_detect(Item, "Fung & Bact ") | str_detect(Item, "Insecticides ")) %>%
  unlist(use.names = FALSE)

specific_pests_model <-list()
coeffs$specific_pests <- tibble()
for (gopi in genes_of_particular_interest) {
  for (item in spec_pests) {
    temp_test_tibble <- test_tibble_log_tpm[[gopi]] %>% 
      rename(est_use_in_2k_radius = all_of(item))
    specific_pests_model[[gopi]][[item]] <- lmer(log_tpm ~ est_use_in_2k_radius + Gut_part + Season +
                                                        (1 | Hive_ID), data = temp_test_tibble)
    
    coeffs$specific_pests <- summary(specific_pests_model[[gopi]][[item]])$coefficients %>%
      as.data.frame() %>%
      rownames_to_column("metric") %>%
      tibble() %>%
      mutate(gene = gopi, 
             Item = item,
             .before = metric) %>%
      mutate(singular = ifelse(isSingular(specific_pests_model[[gopi]][[item]]), TRUE, FALSE)) %>%
      rbind(coeffs$specific_pests)
    
  }
}

#####
# EXTRACT SLOPES
slopes <- list()
for (level in names(coeffs)) {
  slopes[[level]] <- coeffs[[level]] %>%
  filter(metric == "ha_cropland_in_2k_radius" | metric == "est_use_in_2k_radius") %>%
  mutate(raw_p_significant = case_when(`Pr(>|t|)` <= 0.001 ~ "***",
                                       `Pr(>|t|)` <= 0.01 ~ "**",
                                       `Pr(>|t|)` <= 0.05 ~ "*",
                                       `Pr(>|t|)` <= 0.075 ~ ".",
                                       .default = "n.s."
  )) %>%
  mutate(test_name = paste0(gene, "; ", Item), .before = gene)
}

layered_correction_list <- layered_p_adjustments_simple_model(slop = slopes)

lowest_highest <- cropland_and_FAO %>%
  pivot_longer(-Country, names_to = "Item") %>%
  group_by(Item) %>%
  summarise(lowest = min(value),
            highest = max(value))

simple_model_tibble_focused <- layered_correction_list$adjusted_p_values %>%
  filter(p_adjusted <= 0.05) %>%
  inner_join(bind_rows(slopes), ., by = "test_name") %>%
  mutate(adjust_p_significant = case_when(p_adjusted <= 0.001 ~ "***",
                                          p_adjusted <= 0.01 ~ "**",
                                          p_adjusted <= 0.05 ~ "*",
                                          p_adjusted <= 0.075 ~ ".",
                                       .default = "n.s."
  )) %>%
  left_join(., lowest_highest, by = "Item") %>%
  mutate(change_in_range = Estimate * (highest - lowest),
         fold_change_in_range = 10^change_in_range,
         backtrans_estimate = 10^Estimate-1,
         backtrans_error = 10^Estimate - 10^(Estimate - `Std. Error`)) %>%
  # This only changes the axis scale but differently for croopland and pesticides:
  mutate(backtrans_estimate = case_when(Item == "ha_cropland_in_2k_radius" ~ backtrans_estimate * 100, 
                                        .default = backtrans_estimate * 1000),
         backtrans_error = case_when(Item == "ha_cropland_in_2k_radius" ~ backtrans_error * 100,
                                        .default = backtrans_error * 1000),
         unit = case_when(Item == "ha_cropland_in_2k_radius" ~ "(km^2)",
                          .default = "(t)")
         )

slope_plot_simple_model_focused <- simple_model_tibble_focused %>%
  mutate(axis_labels = paste0(Item, " ", unit)) %>%
  mutate(axis_labels = fct_rev(fct_inorder(axis_labels))) %>%
  forest_plot(plot_title = "phosphoadenosine phosphosulfate reductase")

fold_change_in_range_plot <- simple_model_tibble_focused %>%
  mutate(Item = fct_rev(fct_inorder(Item))) %>%
  ggplot(aes(x = Item, y = fold_change_in_range)) +
  geom_col() +
  coord_flip() +
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        axis.text.y=element_blank())

patch_simple_model <- slope_plot_simple_model_focused + fold_change_in_range_plot


simple_model_tibble_all_tests <- layered_correction_list$adjusted_p_values %>%
  inner_join(bind_rows(slopes), ., by = "test_name") %>%
  mutate(axis_labels = fct_rev(fct_inorder(test_name))) %>%
  mutate(backtrans_estimate = 10^Estimate-1,
         backtrans_error = 10^Estimate - 10^(Estimate - `Std. Error`))

slope_plot_simple_model_all_tests <- simple_model_tibble_all_tests %>%
  forest_plot(plot_title = "all tests")
  
#####
# Diagnostics

# plot(model_cropland$`phosphoadenosine phosphosulfate reductase`, which = 1)
# qqnorm(resid(model_cropland$`phosphoadenosine phosphosulfate reductase`))
# 
# plot(total_pest_model$`phosphoadenosine phosphosulfate reductase`$`Pesticides (total)`, which = 1)
# qqnorm(resid(total_pest_model$`phosphoadenosine phosphosulfate reductase`$`Pesticides (total)`))
# 
# plot(pest_groups_model$`phosphoadenosine phosphosulfate reductase`$`Fungicides and Bactericides`, which = 1)
# qqnorm(resid(pest_groups_model$`phosphoadenosine phosphosulfate reductase`$`Fungicides and Bactericides`))
# 
# plot(pest_groups_model$`phosphoadenosine phosphosulfate reductase`$Herbicides, which = 1)
# qqnorm(resid(pest_groups_model$`phosphoadenosine phosphosulfate reductase`$Herbicides))


#####


# # Save files
# system("mkdir -p output/R/gene_content/landuse/")
# ggsave("output/R/gene_content/landuse/simple_model_patch.pdf",
#        patch_simple_model, width = 8, height = 5)
# write_delim(gene_tpm, "output/R/gene_content/landuse/gene_tpm.tsv",
#             delim = "\t")
# write_delim(simple_model_tibble_focused, "output/R/gene_content/landuse/simple_model_tibble_focused.tsv",
#             delim = "\t")
# write_delim(layered_correction_list$adjusted_p_values, "output/R/gene_content/landuse/simple_model_adjusted_p.tsv",
#             delim = "\t")
# 
# write_delim(simple_model_tibble_all_tests, "output/R/gene_content/landuse/simple_model_all_tests.tsv",
#             delim = "\t")
# ggsave("output/R/gene_content/landuse/simple_model_all_tests.pdf",
#        slope_plot_simple_model_all_tests, width = 12, height = 55, limitsize = FALSE)
# 
# for (layer in names(layered_correction_list$subgraph_plots)) {
#   for (gene in names(layered_correction_list$subgraph_plots[[layer]])) {
#     ggsave(paste0("output/R/gene_content/landuse/simple_model_subgraph.", gene, ".", layer,".pdf"),
#            layered_correction_list$subgraph_plots[[layer]][[gene]],
#            width = 20, height = 10)
#   }
# }

