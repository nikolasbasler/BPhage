library(lme4)
library(lmerTest)
library(ggrepel)
library(forcats)
library(tidyverse)
library(patchwork)
library(tidytext)

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

coeffs_tpm_simple <- list()

model_tpm <- list()

#####
# CROPLAND
test_tibble_log_tpm <- list()
coeffs_tpm_simple$cropland <- tibble()
for (goi in genes_of_interest) {
    
  test_tibble_log_tpm[[goi]] <- gene_tpm %>%
    filter(gene == goi) %>%
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
           se_intercept = `Std. Error`) %>%
    select(test_name, intercept, se_intercept) %>%
    full_join(temp_slope_tibble, ., by = "test_name") %>%
    relocate(c(intercept, se_intercept), .after = `Std. Error`)
  
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
  mutate(effect = linear_effect_fun(s = Estimate, h = highest, l = lowest)) %>%
  group_by(gene) %>%
  mutate(effect = ifelse(effect < 1 , 1/effect, effect),
         y_stretching_factor = max(effect) / effect) %>%
  # mutate(y_stretching_factor = max(effect) / effect) %>%
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

color_list <- list(dark = list( `PAPS reductase` = "#ef8f01"),
                   bright = list( `PAPS reductase` = "#ef8f01"))

genes_log_tpm_plots <- list()
for (t_name in unique(sig_tests$test_name)) {
      
  log_tpm_test <- sig_tests %>%
    filter(test_name == t_name)
  
  tested_gene <- log_tpm_test$gene
  tested_item <- log_tpm_test$Item

  genes_log_tpm_plots[[tested_gene]][[tested_item]] <- 
    mixed_model_plot(filt_test_tibble = log_tpm_test,
                     transform_fun = linear_fun,
                     effect_fun = linear_effect_fun,
                     dark_col = color_list$dark[[tested_gene]],
                     bright_col = color_list$bright[[tested_gene]],
                     y_axis_label = "Log rel. gene abund.")
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

# Pest and land use 
facet_order <- c("Pesticides (total)",
              "Insecticides",
              "Herbicides",
              "Fungicides and Bactericides",
              "Plant Growth Regulators",
              spec_pests)
pest_use_tibble <- cropland_and_FAO %>%
  select(Country, Cropland_in_2km_radius, all_of(facet_order)) %>%
  rename(`Cropland in 2km radius` = Cropland_in_2km_radius) %>%
  pivot_longer(-Country, names_to = "parameter") %>%
  mutate(within_ordered_Country = reorder_within(Country, value, parameter),
         parameter = factor(parameter, levels = c("Cropland in 2km radius", facet_order)))

pest_use_plot <- list()
pest_use_plot$country_facet <- pest_use_tibble %>%
  ggplot(aes(x = within_ordered_Country, y = value)) +
  geom_col() +
  geom_text(aes(label = round(value)),
            vjust = -0.3,              
            size = 3) +  
  facet_wrap(~parameter, scales = "free") +
  scale_x_reordered() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(y = "Cropland (ha) or estimated pesticide use (kg) in 2km radius",
       x = "Country")
pest_use_plot$pest_facet <- pest_use_tibble %>%
  ggplot(aes(x = parameter, y = value)) +
  geom_col() +
  geom_text(aes(label = round(value)),
            vjust = -0.3,              
            size = 3) +  
  facet_wrap(~Country, scales = "free") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +
  labs(y = "Cropland (ha) or estimated pesticide use (kg) in 2km radius")

  
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


#####
# SAVE FILES

system("mkdir -p output/R/genes_pathogens_and_landuse/gene_tpm_vs_landuse/single_panels")
system("mkdir -p output/R/genes_pathogens_and_landuse/gene_tpm_vs_landuse/model_diagnostics")

write_delim(bind_rows(coeffs_tpm_simple), "output/R/genes_pathogens_and_landuse/gene_tpm_vs_landuse/gene_tpm_vs_landuse.all_coeffs.tsv",
            delim = "\t")
write_delim(all_slopes, "output/R/genes_pathogens_and_landuse/gene_tpm_vs_landuse/gene_tpm_vs_landuse.all_slopes.tsv",
            delim = "\t")
ggsave("output/R/genes_pathogens_and_landuse/gene_tpm_vs_landuse/gene_tpm_vs_landuse.all_tests.pdf",
       all_tests_forest_plot, width = 12, height = 50, limitsize = FALSE)

write_delim(gene_tpm, "output/R/genes_pathogens_and_landuse/gene_tpm_vs_landuse/gene_tpm.tsv",
            delim = "\t")

ggsave("output/R/genes_pathogens_and_landuse/land_and_pest_use.country_facet.pdf",
       pest_use_plot$country_facet, width = 14, height = 10)
ggsave("output/R/genes_pathogens_and_landuse/land_and_pest_use.pest_facet.pdf",
       pest_use_plot$pest_facet, width = 15, height = 20)
write_delim(pest_use_tibble, "output/R/genes_pathogens_and_landuse/land_and_pest_use.tsv",
            delim = "\t")

for (goi in names(genes_log_tpm_plots)) {
  for (item in names(genes_log_tpm_plots[[goi]])) {
    ggsave(paste0("output/R/genes_pathogens_and_landuse/gene_tpm_vs_landuse/single_panels/", goi, ".", item, ".pdf"),
           genes_log_tpm_plots[[goi]][[item]], height = 3.5, width = 3.5)
    saveRDS(genes_log_tpm_plots[[goi]][[item]],
            paste0("output/R/genes_pathogens_and_landuse/gene_tpm_vs_landuse/single_panels/RDS.", goi, ".", item, ".rds"))
  }
}
ggsave("output/R/genes_pathogens_and_landuse/gene_tpm_vs_landuse/single_panels/common_legend.pdf",
       common_legend, height = 1, width = 6)


for (item in names(model_diagnostics)) {
  for (goi in names(model_diagnostics[[item]])) {
    ggsave(paste0("output/R/genes_pathogens_and_landuse/gene_tpm_vs_landuse/model_diagnostics/model_diagnostics.", item, ".", goi, ".pdf"),
           model_diagnostics[[item]][[goi]], width = 6, height = 6)
  }
}


