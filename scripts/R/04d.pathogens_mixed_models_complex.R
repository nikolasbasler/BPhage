library(lme4)
library(lmerTest)
library(DHARMa)
library(ggrepel)
library(forcats)
library(tidyverse)
library(patchwork)
library(gMCPLite)
library(readxl)

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

pathogen_data <- read_excel("data/GlobalBGOOD_WP1_Tier1_Scien.xlsx", skip = 1) %>%
  rename(BGOOD_sample_code = Sample_ID) %>%
  mutate(CBPV = ifelse(str_detect(CBPV, ">40,00"), "41", CBPV),
         BQCV = as.character(BQCV)) # %>%
# rename(`Cat DWV A` = Cat....24,
#        `Cat DWV B` = Cat....26,
#        `Cat ABPV` = Cat....28,
#        `Cat CBPV` = Cat....30,
#        `Cat BQCV` = Cat....32,
#        `Cat SBV` = Cat....34,
#        `Cat EFB` = Cat....36,
#        `Cat AFB` = Cat....38) %>%
# mutate(across(starts_with("Cat"), ~ na_if(.x, "-")),
#        across(starts_with("Cat"), ~ factor(.x, levels = c("L", "M", "H")))) %>%
# mutate(nosema_spores = ifelse(`N. spores` == "ND", "0", `N. spores`),
#        nosema_spores = ifelse(`N. spores` == '< 25000', "25000" , nosema_spores),
#        nosema_spores = as.numeric(nosema_spores))

pathogens_Cts <- c("DWV A", "DWV B", "ABPV", "CBPV", "BQCV", "SBV", "EFB",
                   "AFB", "N. apis", "N. ceranae")

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


pathogens_of_interest <- c("DWV B", "ABPV", "CBPV", "BQCV", "SBV", "N. ceranae")

coeffs_ct_simple <- list()
countries_left <- list()

#####
# CROPLAND
test_tibble_ct <- list()
model_ct_simple_cropland <- list()
coeffs_ct_simple$cropland <- tibble()
for (poi in pathogens_of_interest) {
  
  test_tibble_ct[[poi]] <- pathogen_ct %>%
    filter(pathogen == poi) %>%
    left_join(., metadata[c("Bee_pool", "Country", "Hive_ID", "Season")], by = "Bee_pool") %>%
    distinct() %>%
    left_join(., cropland_and_FAO, by = "Country") %>%
    mutate(cropland_scaled = scale(ha_cropland_in_2k_radius)[1])
  
  model_ct_simple_cropland[[poi]] <- lmer(Ct ~ ha_cropland_in_2k_radius + Season +
                                            (1 | Country / Hive_ID ), data = test_tibble_ct[[poi]])
  
  has_convergence_issues <- FALSE
  messages <- model_ct_simple_cropland[[poi]]@optinfo$conv$lme4$messages
  if (is.null(messages)) {
    has_convergence_issues <- FALSE
  } else {
    has_convergence_issues <- length(messages[!grepl("singular", messages, ignore.case = TRUE)]) != 0
  }
  if (has_convergence_issues) {
    print(paste0(poi, " - Simple model didn't converge. Removed from the list."))
    model_ct_simple_cropland[[poi]] <- NULL
    
  } else {
    coeffs_ct_simple$cropland <- summary(model_ct_simple_cropland[[poi]])$coefficients %>%
      as.data.frame() %>%
      rownames_to_column("metric") %>%
      tibble() %>%
      mutate(pathogen = poi, .before = metric) %>%
      mutate(Item = metric, .before = metric) %>%
      mutate(singular = ifelse(isSingular(model_ct_simple_cropland[[poi]]), TRUE, FALSE)) %>%
      rbind(coeffs_ct_simple$cropland)
  }
}


##### 
# PESTICIDES:

model_ct_simple_total_pest <- list()
coeffs_ct_simple$total_pest <- tibble()
for (poi in pathogens_of_interest) {
  temp_test_tibble <- test_tibble_ct[[poi]] %>% 
    rename(est_use_in_2k_radius = "Pesticides (total)")
  
  model_ct_simple_total_pest[[poi]][["Pesticides (total)"]] <- lmer(Ct ~ est_use_in_2k_radius + Season +
                                                                      (1 | Country / Hive_ID ), data = temp_test_tibble)
  
  has_convergence_issues <- FALSE
  messages <- model_ct_simple_total_pest[[poi]][["Pesticides (total)"]]@optinfo$conv$lme4$messages
  if (is.null(messages)) {
    has_convergence_issues <- FALSE
  } else {
    has_convergence_issues <- length(messages[!grepl("singular", messages, ignore.case = TRUE)]) != 0
  }
  if (has_convergence_issues) {
    print(paste0(poi, " - Simple model didn't converge. Removed from the list."))
    model_ct_simple_total_pest[[poi]] <- NULL
    
  } else {
    coeffs_ct_simple$total_pest <- summary(model_ct_simple_total_pest[[poi]][["Pesticides (total)"]])$coefficients %>%
      as.data.frame() %>%
      rownames_to_column("metric") %>%
      tibble() %>%
      mutate(pathogen = poi, 
             Item = "Pesticides (total)",
             .before = metric) %>%
      mutate(singular = ifelse(isSingular(model_ct_simple_total_pest[[poi]][["Pesticides (total)"]]), TRUE, FALSE)) %>%
      rbind(coeffs_ct_simple$total_pest)
  }
}


#####
# PEST GROUPS
model_ct_simple_pest_groups <- list()
coeffs_ct_simple$pest_groups <- tibble()
for (poi in pathogens_of_interest) {
  for (item in c("Insecticides", "Herbicides", "Fungicides and Bactericides", "Plant Growth Regulators")) {
    temp_test_tibble <- test_tibble_ct[[poi]] %>% 
      rename(est_use_in_2k_radius = all_of(item))
    model_ct_simple_pest_groups[[poi]][[item]] <- lmer(Ct ~ est_use_in_2k_radius + Season +
                                                         (1 | Country / Hive_ID ), data = temp_test_tibble)
    
    has_convergence_issues <- FALSE
    messages <- model_ct_simple_pest_groups[[poi]][[item]]@optinfo$conv$lme4$messages
    if (is.null(messages)) {
      has_convergence_issues <- FALSE
    } else {
      has_convergence_issues <- length(messages[!grepl("singular", messages, ignore.case = TRUE)]) != 0
    }
    if (has_convergence_issues) {
      print(paste0(poi, " - Simple model didn't converge. Removed from the list."))
      model_ct_simple_pest_groups[[poi]][[item]] <- NULL
      
    } else {
      coeffs_ct_simple$pest_groups <- summary(model_ct_simple_pest_groups[[poi]][[item]])$coefficients %>%
        as.data.frame() %>%
        rownames_to_column("metric") %>%
        tibble() %>%
        mutate(pathogen = poi, 
               Item = item,
               .before = metric) %>%
        mutate(singular = ifelse(isSingular(model_ct_simple_pest_groups[[poi]][[item]]), TRUE, FALSE)) %>%
        rbind(coeffs_ct_simple$pest_groups)
    }
  }
}



#####
# SPECIFIC PESTS

spec_pests <- tibble(Item = colnames(cropland_and_FAO)) %>%
  filter(str_detect(Item, "Herbicides ") | str_detect(Item, "Fung & Bact ") | str_detect(Item, "Insecticides ")) %>%
  unlist(use.names = FALSE)

model_ct_simple_specific_pests <-list()
coeffs_ct_simple$specific_pests <- tibble()
for (poi in pathogens_of_interest) {
  for (item in spec_pests) {
    temp_test_tibble <- test_tibble_ct[[poi]] %>% 
      rename(est_use_in_2k_radius = all_of(item))
    model_ct_simple_specific_pests[[poi]][[item]] <- lmer(Ct ~ est_use_in_2k_radius + Season +
                                                            (1 | Country / Hive_ID ), data = temp_test_tibble)
    
    has_convergence_issues <- FALSE
    messages <- model_ct_simple_specific_pests[[poi]][[item]]@optinfo$conv$lme4$messages
    if (is.null(messages)) {
      has_convergence_issues <- FALSE
    } else {
      has_convergence_issues <- length(messages[!grepl("singular", messages, ignore.case = TRUE)]) != 0
    }
    if (has_convergence_issues) {
      print(paste0(poi, " - Simple model didn't converge. Removed from the list."))
      model_ct_simple_specific_pests[[poi]][[item]] <- NULL
      
    } else {
      coeffs_ct_simple$specific_pests <- summary(model_ct_simple_specific_pests[[poi]][[item]])$coefficients %>%
        as.data.frame() %>%
        rownames_to_column("metric") %>%
        tibble() %>%
        mutate(pathogen = poi, 
               Item = item,
               .before = metric) %>%
        mutate(singular = ifelse(isSingular(model_ct_simple_specific_pests[[poi]][[item]]), TRUE, FALSE)) %>%
        rbind(coeffs_ct_simple$specific_pests)
    }
    
  }
}


#####
# EXTRACT SLOPES
slopes <- list()
for (level in names(coeffs_ct_simple)) {
  slopes[[level]] <- coeffs_ct_simple[[level]] %>%
    filter(metric == "ha_cropland_in_2k_radius" | metric == "est_use_in_2k_radius") %>%
    rename(raw_p_value = `Pr(>|t|)`) %>%
    mutate(raw_p_significant = case_when(raw_p_value <= 0.001 ~ "***",
                                         raw_p_value <= 0.01 ~ "**",
                                         raw_p_value <= 0.05 ~ "*",
                                         raw_p_value <= 0.075 ~ ".",
                                         .default = "n.s."
    )) %>%
    mutate(test_name = paste0(pathogen, "; ", Item), .before = pathogen)
}


##### 
# ADJUST P-VALUES

layered_correction_list <- layered_p_adjustments(slop = slopes, gene_or_pathogen = "pathogen")
layered_correction_list$adjusted_p_values %>% arrange(p_adjusted)

all_slopes <- bind_rows(slopes) %>% 
  left_join(., layered_correction_list$adjusted_p_values, by = "test_name") %>%
  mutate(p_adjust_significant = case_when(p_adjusted <= 0.001 ~ "***",
                                          p_adjusted <= 0.01 ~ "**",
                                          p_adjusted <= 0.05 ~ "*",
                                          p_adjusted <= 0.075 ~ ".",
                                          .default = "n.s."
  ))

# ### TRY OUT:
# all_slopes %>% 
#   mutate(layer = case_when(Item == "ha_cropland_in_2k_radius" ~ "layer_1",
#                            Item == "Pesticides (total)" ~ "layer_2",
#                            Item %in% c("Insecticides", "Herbicides", "Fungicides and Bactericides", "Plant Growth Regulators") ~ "layer_3",
#                            Item %in% spec_pests ~ "layer_4"
#   )) %>%
#   mutate(p_all_adjust = p.adjust(raw_p_value, method = "BH")) %>%
#   mutate(p_all_adjust_significant = case_when(p_all_adjust <= 0.001 ~ "***",
#                                               p_all_adjust <= 0.01 ~ "**",
#                                               p_all_adjust <= 0.05 ~ "*",
#                                               p_all_adjust <= 0.075 ~ ".",
#                                               .default = "n.s."
#   )) %>%
#   group_by(layer) %>%
#   mutate(p_layer_adjust = p.adjust(raw_p_value, method = "BH")) %>%
#   ungroup() %>%
#   mutate(p_layer_adjust_significant = case_when(p_layer_adjust <= 0.001 ~ "***",
#                                                 p_layer_adjust <= 0.01 ~ "**",
#                                                 p_layer_adjust <= 0.05 ~ "*",
#                                                 p_layer_adjust <= 0.075 ~ ".",
#                                                 .default = "n.s."
#   )) %>%
#   filter(p_adjusted <= 0.05 | p_all_adjust <= 0.05 | p_layer_adjust <= 0.05 )  %>%
#   write_delim("output/R/gene_content/landuse/adjust_comparison.pathogens.ct.complex.tsv", delim = "\t")

#####
# MAKE PLOTS

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
  mutate(ct_change_in_range = Estimate * (highest - lowest))

slope_plot_simple_model_focused <- simple_model_tibble_focused %>%
  mutate(axis_labels = test_name) %>%
  mutate(estimate = Estimate,
         error = `Std. Error`) %>%
  forest_plot(plot_title = "pathogens")

fold_change_in_range_plot <- simple_model_tibble_focused %>%
  mutate(Item = fct_rev(fct_inorder(Item))) %>%
  ggplot(aes(x = Item, y = ct_change_in_range)) +
  geom_col() +
  coord_flip() +
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        axis.text.y=element_blank())

patch_simple_model <- slope_plot_simple_model_focused + fold_change_in_range_plot


simple_model_tibble_all_tests <- layered_correction_list$adjusted_p_values %>%
  inner_join(bind_rows(slopes), ., by = "test_name") %>%
  mutate(axis_labels = fct_rev(fct_inorder(test_name))) %>%
  mutate(estimate = Estimate,
         error = `Std. Error`)

slope_plot_simple_model_all_tests <- simple_model_tibble_all_tests %>%
  forest_plot(plot_title = "all tests")

