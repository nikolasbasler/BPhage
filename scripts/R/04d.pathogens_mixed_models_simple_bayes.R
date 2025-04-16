library(lme4)
library(lmerTest)
library(DHARMa)
library(ggrepel)
library(forcats)
library(tidyverse)
library(patchwork)
library(gMCPLite)
library(readxl)
library(brms)

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

# pathogen_ct <- pathogen_data %>% 
#   filter(Years == 2020) %>%
#   select(BGOOD_sample_code, all_of(pathogens_Cts)) %>%
#   left_join(metadata[c("BGOOD_sample_code", "Bee_pool")], ., by = "BGOOD_sample_code") %>%
#   distinct() %>%
#   select(-BGOOD_sample_code) %>%
#   pivot_longer(-Bee_pool, names_to = "pathogen", values_to = "Ct") %>%
#   filter(!is.na(Bee_pool)) %>% 
#   mutate(Ct = str_replace_all(Ct, "negative", "41")) %>%
#   mutate(Ct = as.numeric(Ct)) %>% 
#   filter(!is.na(Ct))

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
  pivot_longer(-Bee_pool, names_to = "pathogen", values_to = "Ct") %>%
  mutate(censoring = case_when(pathogen %in% pathogens & Ct >= 41 ~ "right",
                               pathogen %in% pathogens & Ct < 41 ~ "none",
                               pathogen == "nosema_spores" & Ct <= 25000 & Ct != 0 ~ "left",
                               is.na(Ct) ~ NA,
                               .default = "none"
  ))


pathogens_of_interest <- c("DWV B", "ABPV", "CBPV", "BQCV", "SBV", "N. ceranae", "nosema_spores")

hists <- list()
for (pat in pathogens_of_interest) {
  hists[[pat]] <- pathogen_ct %>%
    filter(pathogen == pat,
           Ct < 41) %>%
    ggplot(aes(x = Ct)) +
    geom_histogram() +
    ggtitle(pat)
}


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
    left_join(., cropland_and_FAO, by = "Country") %>%
    mutate(Ct = if_else(rep(poi == "nosema_spores", n()), as.integer(Ct), Ct))  
  
  # model_ct_simple_cropland[[poi]] <- lmer(Ct ~ ha_cropland_in_2k_radius + Season +
  #                                            ( 1 | Hive_ID ), data = test_tibble_ct[[poi]])
  
  model_ct_censored[[poi]] <- brm(
    formula = bf(Ct | cens(censoring) ~ ha_cropland_in_2k_radius + Season + (1 | Hive_ID )),
    seed = 1,
    data = test_tibble_ct[[poi]],
    prior(normal(20, 5), class = "Intercept", lb = 1, ub = 40),
    family = brmsfamily("gaussian", link = "softplus"),
    # family = gaussian(),
    chains = 4,         # 4 independent MCMC chains
    iter = 10000,        # 2000 iterations per chain (roughly 1000 warmup, 1000 sampling)
    cores = 4           # Utilize 4 CPU cores in parallel
  )
  
  
  summary(model_ct_censored[[poi]])
  
  
  # has_convergence_issues <- FALSE
  # messages <- model_ct_simple_cropland[[poi]]@optinfo$conv$lme4$messages
  # if (is.null(messages)) {
  #   has_convergence_issues <- FALSE
  # } else {
  #   has_convergence_issues <- length(messages[!grepl("singular", messages, ignore.case = TRUE)]) != 0
  # }
  # if (has_convergence_issues) {
  #   print(paste0(poi, " - Simple model didn't converge. Removed from the list."))
  #   model_ct_simple_cropland[[poi]] <- NULL
  #   
  # } else {
  #   coeffs_ct_simple$cropland <- summary(model_ct_simple_cropland[[poi]])$coefficients %>%
  #     as.data.frame() %>%
  #     rownames_to_column("metric") %>%
  #     tibble() %>%
  #     mutate(pathogen = poi, .before = metric) %>%
  #     mutate(Item = metric, .before = metric) %>%
  #     mutate(singular = ifelse(isSingular(model_ct_simple_cropland[[poi]]), TRUE, FALSE)) %>%
  #     rbind(coeffs_ct_simple$cropland)
  # }
}

summary(model_ct_censored$`DWV B`)
summary(model_ct_censored$ABPV)
summary(model_ct_censored$CBPV)
summary(model_ct_censored$BQCV)
summary(model_ct_censored$SBV)
summary(model_ct_censored$`N. ceranae`)
summary(model_ct_censored$nosema_spores)


##### 
# PESTICIDES:

model_ct_simple_total_pest <- list()
coeffs_ct_simple$total_pest <- tibble()
for (poi in pathogens_of_interest) {
  temp_test_tibble <- test_tibble_ct[[poi]] %>% 
    rename(est_use_in_2k_radius = "Pesticides (total)")
  
  model_ct_simple_total_pest[[poi]][["Pesticides (total)"]] <- brm(
    formula = bf(Ct | cens(censoring) ~ est_use_in_2k_radius + Season + (1 | Hive_ID )),
    seed = 1,
    data = temp_test_tibble,
    prior(normal(20, 5), class = "Intercept", lb = 1, ub = 40),
    family = brmsfamily("gaussian", link = "softplus"),
    # family = gaussian(),
    chains = 4,         # 4 independent MCMC chains
    iter = 10000,        # 2000 iterations per chain (roughly 1000 warmup, 1000 sampling)
    cores = 4           # Utilize 4 CPU cores in parallel
  )
}


#####
# PEST GROUPS
model_ct_simple_pest_groups <- list()
coeffs_ct_simple$pest_groups <- tibble()
for (poi in pathogens_of_interest) {
  for (item in c("Insecticides", "Herbicides", "Fungicides and Bactericides", "Plant Growth Regulators")) {
    temp_test_tibble <- test_tibble_ct[[poi]] %>% 
      rename(est_use_in_2k_radius = all_of(item))
    
    model_ct_simple_pest_groups[[poi]][[item]]  <- brm(
      formula = bf(Ct | cens(censoring) ~ est_use_in_2k_radius + Season + (1 | Hive_ID )),
      seed = 1,
      data = temp_test_tibble,
      prior(normal(20, 5), class = "Intercept", lb = 1, ub = 40),
      family = brmsfamily("gaussian", link = "softplus"),
      # family = gaussian(),
      chains = 4,         # 4 independent MCMC chains
      iter = 10000,        # 2000 iterations per chain (roughly 1000 warmup, 1000 sampling)
      cores = 4           # Utilize 4 CPU cores in parallel
    )
    
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

    model_ct_simple_specific_pests[[poi]][[item]] <- brm(
      formula = bf(Ct | cens(censoring) ~ est_use_in_2k_radius + Season + (1 | Hive_ID )),
      seed = 1,
      data = temp_test_tibble,
      prior(normal(20, 5), class = "Intercept", lb = 1, ub = 40),
      family = brmsfamily("gaussian", link = "softplus"),
      # family = gaussian(),
      chains = 4,         # 4 independent MCMC chains
      iter = 10000,        # 2000 iterations per chain (roughly 1000 warmup, 1000 sampling)
      cores = 4           # Utilize 4 CPU cores in parallel
    )
    
  }
}

system("mkdir -p temp.bayes/")
saveRDS(model_ct_censored, "temp.bayes/crop.rds")
saveRDS(model_ct_simple_total_pest, "temp.bayes/total_pest.rds")
saveRDS(model_ct_simple_pest_groups, "temp.bayes/pest_groups.rds")
saveRDS(model_ct_simple_specific_pests, "temp.bayes/spec_pests.rds")




