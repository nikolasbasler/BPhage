library(lme4)
library(lmerTest)
library(DHARMa)
library(ggrepel)
library(forcats)
library(tidyverse)
library(patchwork)
library(gMCPLite)
library(readxl)
library(tidytext)


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
  left_join(FAOSTAT_added_data[c("Country", "cropland_fraction_2k_radius", "Cropland_in_2km_radius")],. , by = "Country") %>%
  distinct() %>%
  select_if(~ !any(is.na(.)))

metadata <- readRDS("output/R/R_variables/metadata.RDS") %>%
  mutate(Hive_ID = as.character(Hive_ID))
classification <- readRDS("output/R/R_variables/classification.RDS")

phage_tpm <- read.csv("output/R/relative_abundance/phage_tpm.csv") %>%
  tibble()

all_hosts <- read.csv("output/R/host_pies/all_hosts.all.csv") %>%
  tibble() %>%
  rename(contig = Virus,
         host_genus = Genus)

hosts_of_interest <- c("Gilliamella", "Snodgrassella", "Bifidobacterium", "Lactobacillus", "Bombilactobacillus",
                       "Frischella", "Bartonella", "Bombella", "Commensalibacter")
host_tpm <- all_hosts %>%
  filter(host_genus %in% hosts_of_interest) %>%
  left_join(., phage_tpm, by = "contig") %>% 
  pivot_longer(-c(contig, host_genus), names_to = "Sample_ID", values_to = "tpm") %>%
  group_by(host_genus, Sample_ID) %>%
  mutate(tpm = sum(tpm)) %>%
  ungroup() %>%
  select(-contig) %>%
  distinct()

coeffs_tpm_simple <- list()


#####
# CROPLAND
test_tibble_log_tpm <- list()
model_tpm_simple_cropland <- list()
coeffs_tpm_simple$cropland <- tibble()
for (hoi in hosts_of_interest) {
  
  test_tibble_log_tpm[[hoi]] <- host_tpm %>%
    filter(host_genus == hoi) %>%
    left_join(., metadata[c("Sample_ID", "Country", "Hive_ID", "Season", "Gut_part")], by = "Sample_ID") %>%
    distinct() %>%
    left_join(., cropland_and_FAO, by = "Country") %>%
    mutate(Gut_part = factor(Gut_part, levels = c("rec", "ile", "mid"))) %>%
    mutate(presence = ifelse(tpm > 0, 1, 0), .before = tpm)

  
  model_tpm_simple_cropland[[hoi]] <- glmer(presence ~ Cropland_in_2km_radius + Gut_part + Season +
                                             (1 | Hive_ID ), data = test_tibble_log_tpm[[hoi]],
                                            family = binomial)
  
  has_convergence_issues <- FALSE
  messages <- model_tpm_simple_cropland[[hoi]]@optinfo$conv$lme4$messages
  if (is.null(messages)) {
    has_convergence_issues <- FALSE
  } else {
    has_convergence_issues <- length(messages[!grepl("singular", messages, ignore.case = TRUE)]) != 0
  }
  if (has_convergence_issues) {
    print(paste0(hoi, " - Simple model didn't converge. Removed from the list."))
    model_tpm_simple_cropland[[hoi]] <- NULL
    
  } else {
    coeffs_tpm_simple$cropland <- summary(model_tpm_simple_cropland[[hoi]])$coefficients %>%
      as.data.frame() %>%
      rownames_to_column("metric") %>%
      tibble() %>%
      mutate(host_genus = hoi, .before = metric) %>%
      mutate(Item = metric, .before = metric) %>%
      mutate(singular = ifelse(isSingular(model_tpm_simple_cropland[[hoi]]), TRUE, FALSE)) %>%
      rbind(coeffs_tpm_simple$cropland)
  }
}



##### 
# PESTICIDES:

model_tpm_simple_total_pest <- list()
coeffs_tpm_simple$total_pest <- tibble()
for (hoi in hosts_of_interest) {
  temp_test_tibble <- test_tibble_log_tpm[[hoi]] %>% 
    rename(est_use_in_2k_radius = "Pesticides (total)")
  model_tpm_simple_total_pest[[hoi]][["Pesticides (total)"]] <- glmer(presence ~ est_use_in_2k_radius + Gut_part + Season +
                                                                       (1 | Hive_ID ), data = temp_test_tibble,
                                                                      family = binomial)
  
  has_convergence_issues <- FALSE
  messages <- model_tpm_simple_total_pest[[hoi]][["Pesticides (total)"]]@optinfo$conv$lme4$messages
  if (is.null(messages)) {
    has_convergence_issues <- FALSE
  } else {
    has_convergence_issues <- length(messages[!grepl("singular", messages, ignore.case = TRUE)]) != 0
  }
  if (has_convergence_issues) {
    print(paste0(hoi, " - Simple model didn't converge. Removed from the list."))
    model_tpm_simple_total_pest[[hoi]] <- NULL
    
  } else {
    coeffs_tpm_simple$total_pest <- summary(model_tpm_simple_total_pest[[hoi]][["Pesticides (total)"]])$coefficients %>%
      as.data.frame() %>%
      rownames_to_column("metric") %>%
      tibble() %>%
      mutate(host_genus = hoi, 
             Item = "Pesticides (total)",
             .before = metric) %>%
      mutate(singular = ifelse(isSingular(model_tpm_simple_total_pest[[hoi]][["Pesticides (total)"]]), TRUE, FALSE)) %>%
      rbind(coeffs_tpm_simple$total_pest)
  }
}

plot(model_tpm_simple_total_pest$Bifidobacterium$`Pesticides (total)`, which = 1)

#####
# PEST GROUPS
model_tpm_simple_pest_groups <- list()
coeffs_tpm_simple$pest_groups <- tibble()
for (hoi in hosts_of_interest) {
  for (item in c("Insecticides", "Herbicides", "Fungicides and Bactericides", "Plant Growth Regulators")) {
    temp_test_tibble <- test_tibble_log_tpm[[hoi]] %>% 
      rename(est_use_in_2k_radius = all_of(item))
    model_tpm_simple_pest_groups[[hoi]][[item]] <- glmer(presence ~ est_use_in_2k_radius + Gut_part + Season +
                                                          (1 | Hive_ID ), data = temp_test_tibble,
                                                         family = binomial)
    
    has_convergence_issues <- FALSE
    messages <- model_tpm_simple_pest_groups[[hoi]][[item]]@optinfo$conv$lme4$messages
    if (is.null(messages)) {
      has_convergence_issues <- FALSE
    } else {
      has_convergence_issues <- length(messages[!grepl("singular", messages, ignore.case = TRUE)]) != 0
    }
    if (has_convergence_issues) {
      print(paste0(hoi, " - Simple model didn't converge. Removed from the list."))
      model_tpm_simple_pest_groups[[hoi]][[item]] <- NULL
      
    } else {
      coeffs_tpm_simple$pest_groups <- summary(model_tpm_simple_pest_groups[[hoi]][[item]])$coefficients %>%
        as.data.frame() %>%
        rownames_to_column("metric") %>%
        tibble() %>%
        mutate(host_genus = hoi, 
               Item = item,
               .before = metric) %>%
        mutate(singular = ifelse(isSingular(model_tpm_simple_pest_groups[[hoi]][[item]]), TRUE, FALSE)) %>%
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
for (hoi in hosts_of_interest) {
  for (item in spec_pests) {
    temp_test_tibble <- test_tibble_log_tpm[[hoi]] %>% 
      rename(est_use_in_2k_radius = all_of(item))
    model_tpm_simple_specific_pests[[hoi]][[item]] <- glmer(presence ~ est_use_in_2k_radius + Gut_part + Season +
                                                             (1 | Hive_ID ), data = temp_test_tibble,
                                                            family = binomial)
    
    has_convergence_issues <- FALSE
    messages <- model_tpm_simple_specific_pests[[hoi]][[item]]@optinfo$conv$lme4$messages
    if (is.null(messages)) {
      has_convergence_issues <- FALSE
    } else {
      has_convergence_issues <- length(messages[!grepl("singular", messages, ignore.case = TRUE)]) != 0
    }
    if (has_convergence_issues) {
      print(paste0(hoi, " - Simple model didn't converge. Removed from the list."))
      model_tpm_simple_specific_pests[[hoi]][[item]] <- NULL
      
    } else {
      coeffs_tpm_simple$specific_pests <- summary(model_tpm_simple_specific_pests[[hoi]][[item]])$coefficients %>%
        as.data.frame() %>%
        rownames_to_column("metric") %>%
        tibble() %>%
        mutate(host_genus = hoi, 
               Item = item,
               .before = metric) %>%
        mutate(singular = ifelse(isSingular(model_tpm_simple_specific_pests[[hoi]][[item]]), TRUE, FALSE)) %>%
        rbind(coeffs_tpm_simple$specific_pests)
    }
    
  }
}




#####
# EXTRACT SLOPES
slopes <- list()
for (level in names(coeffs_tpm_simple)) {
  
  if (nrow(coeffs_tpm_simple[[level]]) == 0) {
    next
  }
  
  temp_slope_tibble <- coeffs_tpm_simple[[level]] %>%
    filter(metric %in% c("Cropland_in_2km_radius", "est_use_in_2k_radius")) %>%
    rename(raw_p_value = `Pr(>|z|)`) %>%
    mutate(raw_p_significant = case_when(raw_p_value <= 0.001 ~ "***",
                                         raw_p_value <= 0.01 ~ "**",
                                         raw_p_value <= 0.05 ~ "*",
                                         raw_p_value <= 0.075 ~ ".",
                                         .default = "n.s."
    )) %>%
    mutate(test_name = paste0(host_genus, "; ", Item), .before = host_genus)
  
  slopes[[level]] <- coeffs_tpm_simple[[level]] %>%
    filter(metric == "(Intercept)") %>%
    mutate(test_name = paste0(host_genus, "; ", temp_slope_tibble$Item), .before = host_genus) %>%
    rename(intercept = Estimate,
           sd_intercept = `Std. Error`) %>%
    select(test_name, intercept, sd_intercept) %>%
    full_join(temp_slope_tibble, ., by = "test_name") %>%
    relocate(c(intercept, sd_intercept), .after = `Std. Error`)
  
}

##### 
# ADJUST P-VALUES

# layered_correction_list <- layered_p_adjustments(slop = slopes)

# all_slopes <- bind_rows(slopes) %>% 
#   left_join(., layered_correction_list$adjusted_p_values, by = "test_name") %>%
#   mutate(p_adjust_significant = case_when(p_adjusted <= 0.001 ~ "***",
#                                           p_adjusted <= 0.01 ~ "**",
#                                           p_adjusted <= 0.05 ~ "*",
#                                           p_adjusted <= 0.075 ~ ".",
#                                           .default = "n.s."
#   ))

all_slopes <- bind_rows(slopes) %>%
  mutate(p_adjusted = p.adjust(raw_p_value, method = "BH")) %>%
  mutate(p_adjust_significant = case_when(p_adjusted <= 0.001 ~ "***",
                                          p_adjusted <= 0.01 ~ "**",
                                          p_adjusted <= 0.05 ~ "*",
                                          p_adjusted <= 0.075 ~ ".",
                                          .default = "n.s.")
  )
