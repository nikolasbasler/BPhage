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
  left_join(FAOSTAT_added_data[c("Country", "cropland_fraction_2k_radius", "Cropland_in_2km_radius")],. , by = "Country") %>%
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

countries_with_presence <- list()
samples_with_presence <- list()


#####
# PREVALENCE PLOT

test_tibble_logit <- list()
for (poi in pathogens_Cts) {
  # poi = "N. ceranae"
  
  test_tibble_logit[[poi]] <- pathogen_ct %>%
    filter(pathogen == poi) %>%
    left_join(., metadata[c("Bee_pool", "Country", "Hive_ID", "Season")], by = "Bee_pool") %>%
    distinct() %>%
    left_join(., cropland_and_FAO, by = "Country") %>%
    mutate(presence = ifelse(Ct < 40, 1, 0), .before = Ct)
}


countries_with_presence$cropland <- list()
samples_with_presence$cropland <- list()
for (logit_test in pathogens_Cts) {
  countries_with_presence$cropland <- test_tibble_logit[[logit_test]] %>%
    filter(presence > 0 ) %>%
    distinct(Country) %>%
    reframe(test = logit_test, countries_with_presence = n()) %>%
    rbind(countries_with_presence$cropland, .)
  samples_with_presence$cropland <- test_tibble_logit[[logit_test]] %>%
    filter(presence > 0 ) %>%
    reframe(test = logit_test, samples_with_presence = n()) %>%
    rbind(samples_with_presence$cropland, .)
}

prevalence_plot <- bind_rows(test_tibble_logit) %>%
  select(pathogen, Bee_pool, Country, presence) %>%
  ggplot(aes(x = reorder(pathogen, -presence, FUN = mean), y = presence)) +
  geom_bar(stat = "summary", fun = mean) +
  geom_text(
    stat = "summary",
    fun = mean,
    aes(label = round(after_stat(y), 2)),
    vjust = -0.5,
    # size = 3  # smaller text (default is ~5)
  ) +  
  labs(x = "Pathogen", y = "Prevalence") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +  # add 10% space at the top
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
prevalence_plot_facet_countries <- prevalence_plot + facet_wrap(~Country)

prevalence_plot_facet_pathogens <- bind_rows(test_tibble_logit) %>%
  select(pathogen, Bee_pool, Country, presence) %>%
  ggplot(aes(x = Country, y = presence)) +
  geom_bar(stat = "summary", fun = mean) +
  geom_text(
    stat = "summary",
    fun = mean,
    aes(label = round(after_stat(y), 2)),
    vjust = -0.5,
    # size = 3  # smaller text (default is ~5)
  ) +  
  labs(x = "Pathogen", y = "Prevalence") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +  # add 10% space at the top
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_wrap(~pathogen)

countries_with_presence 
samples_with_presence
prevalence_plot
prevalence_plot_facet_pathogens

###### 
# CROPLAND

# pathogens_of_interest <- c("DWV A", "ABPV", "CBPV", "N. ceranae")
pathogens_of_interest <- c("ABPV", "N. ceranae", "CBPV")

coeffs_logit <- list()
model_logit_cropland <- list()
for (poi in pathogens_of_interest) {

  model_logit_cropland[[poi]] <- glmer(
    presence ~ Cropland_in_2km_radius + Season + ( 1 | Hive_ID ),
    data = test_tibble_logit[[poi]],
    family = binomial)
  
  test_tibble_logit[[poi]] %>% count(Country) %>% arrange(desc(n))
  test_tibble_logit[[poi]] %>% count(Hive_ID) %>% arrange(desc(n))
  test_tibble_logit[[poi]] %>% count(Season) %>% arrange(desc(n))
  # summary(model_logit_cropland[[poi]])

  has_convergence_issues <- FALSE
  messages <- model_logit_cropland[[poi]]@optinfo$conv$lme4$messages
  if (!is.null(messages) & any(str_detect(messages, "failed to converge"))) {
    has_convergence_issues <- TRUE
  }

  if (has_convergence_issues) {
    print(paste0(poi, " - Model didn't converge. Removed from the list."))
    model_logit_cropland[[poi]] <- NULL
  } else {
    coeffs_logit$cropland <- summary(model_logit_cropland[[poi]])$coefficients %>%
      as.data.frame() %>%
      rownames_to_column("metric") %>%
      tibble() %>%
      mutate(pathogen = poi, .before = metric) %>%
      mutate(Item = metric, .before = metric) %>%
      mutate(singular = ifelse(isSingular(model_logit_cropland[[poi]]), TRUE, FALSE)) %>%
      rbind(coeffs_logit$cropland)
  }
}
# 
# test_tibble_joint_logit <- pathogen_ct %>%
#   mutate(presence = ifelse(Ct < 41, 1, 0)) %>%
#   filter(pathogen %in% pathogens_of_interest) %>%
#   select(-Ct) %>%
#   pivot_wider(names_from = pathogen, values_from = presence) %>%
#   filter(if_all(everything(), ~ !is.na(.))) %>%
#   left_join(., metadata[c("Bee_pool", "Country", "Hive_ID", "Season")], by = "Bee_pool") %>%
#   distinct() %>%
#   left_join(., cropland_and_FAO, by = "Country") %>%
#   rename(n_cera = `N. ceranae`)
# 
# test_tibble_joint_logit %>% count(Season) %>% arrange(desc(n))
# 
# model_combinded_logit_cropland <- glmer(
#   ABPV ~ Cropland_in_2km_radius + Season + ( 1 | Hive_ID ),
#   data = test_tibble_joint_logit,
#   family = binomial)
# 
# test_tibble_combinded_logit <- bind_rows(test_tibble_logit) %>%
#   filter(pathogen %in% pathogens_of_interest) %>%
#   group_by(Bee_pool) %>%
#   mutate(num_pathogens = sum(presence), .before = presence) %>%
#   ungroup() %>%
#   select(-c(Ct, presence, pathogen)) %>%
#   distinct()
# 
# test_tibble_combinded_logit %>% count(Country)
# test_tibble_combinded_logit %>% count(Season) %>% arrange(n)

# model_combinded_logit_cropland <- glmer(
#   cbind(num_pathogens, 3 - num_pathogens) ~ Cropland_in_2km_radius + Season + ( 1 | Hive_ID ),
#   data = test_tibble_combinded_logit,
#   family = binomial)
# 
# summary(model_combinded_logit_cropland)


##### 
# PESTICIDES:

model_logit_total_pest <- list()
# coeffs_logit$total_pest <- tibble()
coeffs_logit$total_pest <- coeffs_logit$cropland %>% filter(FALSE)
for (poi in pathogens_of_interest) {
  temp_test_tibble <- test_tibble_logit[[poi]] %>%
    rename(est_use_in_2k_radius = `Pesticides (total)`)
  
  model_logit_total_pest[[poi]][["Pesticides (total)"]] <- glmer(presence ~ est_use_in_2k_radius + Season +
                                                                   ( 1 | Country / Hive_ID ), data = temp_test_tibble,
                                                                 family = binomial)
  # summary( model_logit_total_pest[[poi]][["Pesticides (total)"]])
  
  has_convergence_issues <- FALSE
  messages <- model_logit_total_pest[[poi]][["Pesticides (total)"]]@optinfo$conv$lme4$messages
  if (!is.null(messages) & any(str_detect(messages, "failed to converge"))) {
    has_convergence_issues <- TRUE
  }
  
  if (has_convergence_issues) {
    print(paste0(poi, " - Model didn't converge. Removed from the list."))
    model_logit_total_pest[[poi]] <- NULL
  } else {
    coeffs_logit$total_pest <- summary(model_logit_total_pest[[poi]][["Pesticides (total)"]])$coefficients %>%
      as.data.frame() %>%
      rownames_to_column("metric") %>%
      tibble() %>%
      mutate(pathogen = poi, 
             Item = "Pesticides (total)",
             .before = metric) %>%
      mutate(singular = ifelse(isSingular(model_logit_total_pest[[poi]][["Pesticides (total)"]]), TRUE, FALSE)) %>%
      rbind(coeffs_logit$total_pest)
  }
}


#####
# PEST GROUPS
model_logit_pest_groups <- list()
coeffs_logit$pest_groups <- tibble()
for (poi in pathogens_of_interest) {
  for (item in c("Insecticides", "Herbicides", "Fungicides and Bactericides", "Plant Growth Regulators")) {
    temp_test_tibble <- test_tibble_logit[[poi]] %>% 
      rename(est_use_in_2k_radius = all_of(item))
    
    model_logit_pest_groups[[poi]][[item]] <- glmer(presence ~ est_use_in_2k_radius + Season +
                                                      ( 1 | Country / Hive_ID ), data = temp_test_tibble,
                                                    family = binomial)
    
    has_convergence_issues <- FALSE
    messages <- model_logit_pest_groups[[poi]][[item]]@optinfo$conv$lme4$messages
    if (!is.null(messages) & any(str_detect(messages, "failed to converge"))) {
      has_convergence_issues <- TRUE
    }
    
    if (has_convergence_issues) {
      print(paste0(poi, "; ", item, " - Model didn't converge. Removed from the list."))
      model_logit_pest_groups[[poi]][[item]] <- NULL
    } else {
      coeffs_logit$pest_groups <- summary(model_logit_pest_groups[[poi]][[item]])$coefficients %>%
        as.data.frame() %>%
        rownames_to_column("metric") %>%
        tibble() %>%
        mutate(pathogen = poi, 
               Item = item,
               .before = metric) %>%
        mutate(singular = ifelse(isSingular(model_logit_pest_groups[[poi]][[item]]), TRUE, FALSE)) %>%
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
for (poi in pathogens_of_interest) {
  for (item in spec_pests) {
    temp_test_tibble <- test_tibble_logit[[poi]] %>% 
      rename(est_use_in_2k_radius = all_of(item))
    
    model_logit_specific_pests[[poi]][[item]] <- glmer(presence ~ est_use_in_2k_radius + Season +
                                                         ( 1 | Country / Hive_ID ), data = temp_test_tibble,
                                                       family = binomial)
    
    has_convergence_issues <- FALSE
    messages <- model_logit_specific_pests[[poi]][[item]]@optinfo$conv$lme4$messages
    if (!is.null(messages) & any(str_detect(messages, "failed to converge"))) {
      has_convergence_issues <- TRUE
    }
    if (has_convergence_issues) {
      print(paste0(poi, "; ", item, " - Model didn't converge. Removed from the list."))
      model_logit_pest_groups[[poi]][[item]] <- NULL
    } else {
      coeffs_logit$specific_pests <- summary(model_logit_specific_pests[[poi]][[item]])$coefficients %>%
        as.data.frame() %>%
        rownames_to_column("metric") %>%
        tibble() %>%
        mutate(pathogen = poi, 
               Item = item,
               .before = metric) %>%
        mutate(singular = ifelse(isSingular(model_logit_specific_pests[[poi]][[item]]), TRUE, FALSE)) %>%
        rbind(coeffs_logit$specific_pests)
    }
  }
}


#####
# EXTRACT SLOPES
slopes <- list()
for (level in names(coeffs_logit)) {
  slopes[[level]] <- coeffs_logit[[level]] %>%
    filter(metric == "Cropland_in_2km_radius" | metric == "est_use_in_2k_radius") %>%
    rename(raw_p_value = `Pr(>|z|)`) %>%
    mutate(raw_p_significant = case_when(raw_p_value <= 0.001 ~ "***",
                                         raw_p_value <= 0.01 ~ "**",
                                         raw_p_value <= 0.05 ~ "*",
                                         raw_p_value <= 0.075 ~ ".",
                                         .default = "n.s."
    )) %>%
    mutate(test_name = paste0(pathogen, "; ", Item), .before = pathogen)
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
# SAVE FILES
system("mkdir -p output/R/gene_content/landuse/pathogen_logit_model")
write_delim(all_slopes, "output/R/gene_content/landuse/pathogen_logit_model/pathogens_logit_model_tibble.tsv",
            delim = "\t")
ggsave("output/R/gene_content/landuse/pathogen_logit_model/pathogens_prevalence_total.pdf",
       prevalence_plot, width = 5, height = 5)
ggsave("output/R/gene_content/landuse/pathogen_logit_model/pathogens_prevalence_country_facet.pdf",
       prevalence_plot_facet_countries, width = 10, height = 10)
ggsave("output/R/gene_content/landuse/pathogen_logit_model/pathogens_prevalence_gene_facet.pdf",
       prevalence_plot_facet_pathogens, width = 10, height = 8)


