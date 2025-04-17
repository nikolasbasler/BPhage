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

pathogens <- c("DWV A", "DWV B", "ABPV", "CBPV", "BQCV", "SBV", "EFB",
               "AFB", "N. apis", "N. ceranae", "N. spores")

pathogen_ct <- pathogen_data %>% 
  filter(Years == 2020) %>%
  select(BGOOD_sample_code, all_of(pathogens_Cts), "N. spores") %>%
  mutate(nosema_spores = ifelse(`N. spores` == "ND", "0", `N. spores`),
         nosema_spores = ifelse(`N. spores` == '< 25000', "25000" , nosema_spores),
         nosema_spores = as.numeric(nosema_spores)) %>%
  left_join(metadata[c("BGOOD_sample_code", "Bee_pool")], ., by = "BGOOD_sample_code") %>%
  distinct() %>% 
  tibble() %>%
  filter(!is.na(Bee_pool)) %>%
  mutate(across(all_of(pathogens_Cts), ~ str_replace_all(.x, "negative", "41")),
         across(all_of(pathogens_Cts), ~ as.numeric(.x)),
         # across(all_of(pathogens_Cts), ~ 2^-.x),
         # across(all_of(pathogens_Cts), ~ ifelse(.x < 2^-40, 0, .x))
         ) %>%
  select(-c(BGOOD_sample_code, `N. spores`)) %>%
  pivot_longer(-Bee_pool, names_to = "pathogen", values_to = "Ct")


# pathogens_of_interest <- c("DWV B", "ABPV", "CBPV", "BQCV", "SBV", "N. ceranae", "nosema_spores")
pathogens_of_interest <- c("DWV B", "BQCV", "SBV")

coeffs_ct_simple <- list()
countries_left <- list()

#####
# CROPLAND
test_tibble_ct <- list()
model_ct_simple_cropland <- list()
coeffs_ct_simple$cropland <- tibble()
model_ct_censored <- list()
for (poi in pathogens_of_interest) {
  
  test_tibble_ct[[poi]] <- pathogen_ct %>%
    filter(pathogen == poi) %>%
    left_join(., metadata[c("Bee_pool", "Country", "Hive_ID", "Season")], by = "Bee_pool") %>%
    distinct() %>%
    left_join(., cropland_and_FAO, by = "Country")
    # filter(Ct < 40)
    # mutate(Ct = if_else(rep(poi == "nosema_spores", n()), as.integer(Ct), Ct))  

  model_ct_simple_cropland[[poi]] <- lmer(Ct ~ ha_cropland_in_2k_radius + Season +
                                             ( 1 | Hive_ID ), data = test_tibble_ct[[poi]])
  
  

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

summary(model_ct_simple_cropland$`DWV B`)
summary(model_ct_simple_cropland$BQCV)
summary(model_ct_simple_cropland$SBV)

# plot(model_ct_simple_cropland$`DWV B`, which = 1)
# qqnorm(resid(model_ct_simple_cropland$`DWV B`))
# 
# plot(model_ct_simple_cropland$BQCV, which = 1)
# qqnorm(resid(model_ct_simple_cropland$BQCV))
# 
# plot(model_ct_simple_cropland$SBV, which = 1)
# qqnorm(resid(model_ct_simple_cropland$SBV))

# 
# test_tibble_combinded_ct <- pathogen_ct %>%
#   filter(pathogen %in% pathogens_of_interest,
#          Ct < 41) %>%
#   mutate(rel_load = 2^-Ct) %>%
#   group_by(Bee_pool) %>%
#   mutate(summed_rel_load = sum(rel_load),
#          summed_back_ct = -log2(summed_rel_load)) %>%
#   ungroup() %>%
#   select(-c(pathogen, rel_load, summed_rel_load, Ct)) %>%
#   distinct() %>%
#   left_join(., metadata[c("Bee_pool", "Country", "Hive_ID", "Season")], by = "Bee_pool") %>%
#   distinct() %>%
#   left_join(., cropland_and_FAO, by = "Country") 
#   # mutate(Ct = if_else(rep(poi == "nosema_spores", n()), as.integer(Ct), Ct)) 
# 
# model_combinded_ct_simple_cropland <- lmer(summed_back_ct ~ ha_cropland_in_2k_radius + Season +
#                                              ( 1 | Hive_ID ), data = test_tibble_combinded_ct)
# summary(model_combinded_ct_simple_cropland)
#   
# plot(model_combinded_ct_simple_cropland, which = 1)
# qqnorm(resid(model_combinded_ct_simple_cropland))


##### 
# PESTICIDES:

model_ct_simple_total_pest <- list()
coeffs_ct_simple$total_pest <- tibble()
for (poi in pathogens_of_interest) {
  temp_test_tibble <- test_tibble_ct[[poi]] %>% 
    rename(est_use_in_2k_radius = "Pesticides (total)")
  
  model_ct_simple_total_pest[[poi]][["Pesticides (total)"]] <- lmer(Ct ~ est_use_in_2k_radius + Season +
                                                                        ( 1 | Hive_ID ), data = temp_test_tibble)
  
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
                                                           ( 1 | Hive_ID ), data = temp_test_tibble)
    
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
                                                              ( 1 | Hive_ID  ), data = temp_test_tibble)
    
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
  
  temp_slope_tibble <- coeffs_ct_simple[[level]] %>%
    filter(metric == "ha_cropland_in_2k_radius" | metric == "est_use_in_2k_radius") %>%
    rename(raw_p_value = `Pr(>|t|)`) %>%
    mutate(raw_p_significant = case_when(raw_p_value <= 0.001 ~ "***",
                                         raw_p_value <= 0.01 ~ "**",
                                         raw_p_value <= 0.05 ~ "*",
                                         raw_p_value <= 0.075 ~ ".",
                                         .default = "n.s."
    )) %>%
    mutate(test_name = paste0(pathogen, "; ", Item), .before = pathogen)
  
  slopes[[level]] <- coeffs_ct_simple[[level]] %>%
    filter(metric == "(Intercept)") %>%
    mutate(test_name = paste0(pathogen, "; ", temp_slope_tibble$Item), .before = pathogen) %>%
    rename(intercept = Estimate,
           sd_intercept = `Std. Error`) %>%
    select(test_name, intercept, sd_intercept) %>%
    full_join(temp_slope_tibble, ., by = "test_name") %>%
    relocate(c(intercept, sd_intercept), .after = `Std. Error`)


}

all_slopes <- bind_rows(slopes) %>%
  mutate(p_adjusted = p.adjust(raw_p_value, method = "BH")) %>%
  mutate(p_adjust_significant = case_when(p_adjusted <= 0.001 ~ "***",
                                          p_adjusted <= 0.01 ~ "**",
                                          p_adjusted <= 0.05 ~ "*",
                                          p_adjusted <= 0.075 ~ ".",
                                          .default = "n.s.")
  ) %>%
  mutate(layer = ifelse(Item %in% spec_pests, "spec_pests", "upper_layer"),
         layer = factor(layer, levels = c("spec_pests", "upper_layer")))

#####
# MAKE PLOTS

lowest_highest <- cropland_and_FAO %>%
  pivot_longer(-Country, names_to = "Item") %>%
  group_by(Item) %>%
  summarise(lowest = min(value),
            highest = max(value))

sig_tests <- all_slopes %>%
  filter(p_adjusted < 0.05) %>%
  left_join(., lowest_highest, by = "Item") %>%
  mutate(effect = ct_effect_fun(s = Estimate, h = highest, l = lowest)) %>%
  group_by(pathogen) %>%
  mutate(y_stretching_factor = max(abs(effect)) / abs(effect),
         y_stretching_factor = ifelse(effect < 0, 1 / y_stretching_factor, y_stretching_factor)) %>%
  # mutate(y_stretching_factor = max(effect) / effect) %>%
  # mutate(which_y_end_to_stretch = ifelse(pathogen == "DWV B", "lower_end", "upper_end")) %>%
  mutate(which_y_end_to_stretch = if_else(
    Estimate[y_stretching_factor == 1] > 0,
    "upper_end",
    "lower_end")) %>%
  ungroup() %>%
  mutate(y_stretching_factor = case_when(
    which_y_end_to_stretch == "upper_end" & Estimate < 0 ~ 1/y_stretching_factor,
    which_y_end_to_stretch == "lower_end" & Estimate > 0 ~ 1/y_stretching_factor,
    .default = y_stretching_factor)) %>%
  select(-effect)

color_list <- list(dark = list(BQCV = "black", `DWV B` = "#ef8f01"),
                   bright = list(BQCV = "black", `DWV B` = "#8B4513"))


pathogens_ct_plots <- list()
for (t_name in unique(sig_tests$test_name)) {
  
  ct_test <- sig_tests %>%
    filter(test_name == t_name)
  
  tested_pathogen <- ct_test$pathogen
  tested_item <- ct_test$Item
  
  pathogens_ct_plots[[tested_pathogen]][[tested_item]] <- 
    mixed_model_plot(filt_test_tibble = ct_test,
                     transform_fun = linear_fun,
                     effect_fun = ct_effect_fun,
                     dark_col = color_list$dark[[tested_pathogen]],
                     bright_col = color_list$bright[[tested_pathogen]],
                     y_axis_label = "Ct")
}


common_legend <- legend_factory(title = "Pathogen", 
                                items = names(color_list$dark),
                                colors = unlist(color_list$dark),
                                position = "bottom")

wrap_of_wraps <- wrap_plots(
  wrap_plots(pathogens_ct_plots$`DWV B`, nrow = 1, axes = "collect"),
    wrap_plots(pathogens_ct_plots$BQCV, nrow = 1, axes = "collect"),
    common_legend,
  nrow = 3, heights = c(rep(4, 2), 1)
)


# simple_model_tibble_focused <- all_slopes %>%
#   filter(p_adjusted <= 0.05) %>%
#   mutate(adjust_p_significant = case_when(p_adjusted <= 0.001 ~ "***",
#                                           p_adjusted <= 0.01 ~ "**",
#                                           p_adjusted <= 0.05 ~ "*",
#                                           p_adjusted <= 0.075 ~ ".",
#                                           .default = "n.s."
#   )) %>%
#   left_join(., lowest_highest, by = "Item") %>%
#   mutate(ct_change_in_range = Estimate * (highest - lowest))
# 
# slope_plot_simple_model_focused <- simple_model_tibble_focused %>%
#   arrange(test_name) %>%
#   mutate(axis_labels = test_name,
#          axis_labels = fct_inorder(axis_labels)) %>%
#   mutate(estimate = Estimate,
#          error = `Std. Error`) %>%
#   forest_plot(plot_title = "pathogens")
# 
# fold_change_in_range_plot <- simple_model_tibble_focused %>%
#   arrange(test_name) %>%
#   mutate(axis_labels = test_name,
#          axis_labels = fct_inorder(axis_labels)) %>%
#   ggplot(aes(x = axis_labels, y = ct_change_in_range)) +
#   geom_col() +
#   coord_flip() +
#   theme_minimal() +
#   theme(axis.title.y = element_blank(),
#         axis.text.y=element_blank())
# 
# patch_simple_model <- slope_plot_simple_model_focused + fold_change_in_range_plot

simple_model_tibble_all_tests <- all_slopes %>%
  mutate(axis_labels = fct_rev(fct_inorder(test_name))) %>%
  mutate(estimate = 10^Estimate-1,
         error = 10^Estimate - 10^(Estimate - `Std. Error`))

slope_plot_simple_model_all_tests <- simple_model_tibble_all_tests %>%
  forest_plot(plot_title = "all tests")


#####
# DIAGNOSTICS

# plot(model_ct_simple_cropland$`DWV B`, which = 1)
# qqnorm(resid(model_ct_simple_cropland$`DWV B`))
# 
# plot(model_ct_simple_specific_pests$`DWV B`$`Herbicides – Dinitroanilines` , which = 1)
# qqnorm(resid(model_ct_simple_specific_pests$`DWV B`$`Herbicides – Dinitroanilines`))
# 
# plot(model_ct_simple_pest_groups$BQCV$Insecticides, which = 1)
# qqnorm(resid(model_ct_simple_pest_groups$BQCV$Insecticides))
# 
# plot(model_ct_simple_specific_pests$BQCV$`Insecticides - nes`, which = 1)
# qqnorm(resid(model_ct_simple_specific_pests$BQCV$`Insecticides - nes`))


#####
# SAVE FILES

system("mkdir -p output/R/gene_content/landuse/pathogen_simple_model")

ggsave("output/R/gene_content/landuse/pathogen_simple_model/ct_mixed_model_wrap.pdf",
       wrap_of_wraps, width = 6, height = 6)


write_delim(pathogen_ct, "output/R/gene_content/landuse/pathogen_ct.tsv",
            delim = "\t")

write_delim(all_slopes, "output/R/gene_content/landuse/pathogen_simple_model/pathogen_simple_model_all_slopes.tsv",
            delim = "\t")

write_delim(simple_model_tibble_all_tests, "output/R/gene_content/landuse/pathogen_simple_model/pathogen_simple_model_all_tests.tsv",
            delim = "\t")
ggsave("output/R/gene_content/landuse/pathogen_simple_model/pathogen_simple_model_all_tests.pdf",
       slope_plot_simple_model_all_tests, width = 12, height = 30, limitsize = FALSE)

