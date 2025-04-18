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


pathogens_of_interest <- "nosema_spores"

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
  
    model_ct_simple_cropland[[poi]] <- glmer.nb(Ct ~ Cropland_in_2km_radius + Season +
                                              ( 1 | Hive_ID ), data = test_tibble_ct[[poi]])
    
    summary(model_ct_simple_cropland[[poi]])
    plot(model_ct_simple_cropland$nosema_spores, which = 1)
    qqnorm(resid(model_ct_simple_cropland$nosema_spores))
    
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


### APPROACH ABANDONED