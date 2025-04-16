# library(lme4)
# library(lmerTest)
# library(DHARMa)
# library(ggrepel)
# library(forcats)
library(tidyverse)
library(patchwork)
# library(gMCPLite)
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



coeffs_ct_simple <- list()
countries_left <- list()

spec_pests <- tibble(Item = colnames(cropland_and_FAO)) %>%
  filter(str_detect(Item, "Herbicides ") | str_detect(Item, "Fung & Bact ") | str_detect(Item, "Insecticides ")) %>%
  unlist(use.names = FALSE)

#####
# CROPLAND

# pathogens_of_interest <- c("DWV B", "ABPV", "CBPV", "BQCV", "SBV", "N. ceranae", "nosema_spores")
pathogens_of_interest <- c("DWV B", "ABPV", "CBPV", "BQCV", "SBV", "N. ceranae")

test_tibble_ct <- list()
joint_model <- list()
for (poi in pathogens_of_interest) {
# for (poi in "DWV B") {
  
  poi = "N. ceranae"
  
  test_tibble_ct[[poi]] <- pathogen_ct %>%
    filter(pathogen == poi) %>%
    left_join(., metadata[c("Bee_pool", "Country", "Hive_ID", "Season")], by = "Bee_pool") %>%
    distinct() %>%
    left_join(., cropland_and_FAO, by = "Country") %>%
    mutate(Ct = if_else(rep(poi == "nosema_spores", n()), as.integer(Ct), Ct)) %>%
    rename(total_pest = `Pesticides (total)`,
           fung_and_bact = `Fungicides and Bactericides`,
           plant_G = `Plant Growth Regulators`,
           ins_org = `Insecticides – Organo-phosphates`,
           ins_car = `Insecticides – Carbamates`,
           ins_pyr = `Insecticides – Pyrethroids`,
           ins_nes = `Insecticides - nes`,
           
           herb_phe = `Herbicides – Phenoxy hormone products`,
           herb_tri = `Herbicides – Triazines`,
           herb_ami = `Herbicides – Amides`,
           herb_car = `Herbicides – Carbamates`,
           herb_din = `Herbicides – Dinitroanilines`,
           herb_ure = `Herbicides – Urea derivates`,
           herb_nes = `Herbicides - nes`,
           
           fung_ino = `Fung & Bact – Inorganics`,
           fung_dit = `Fung & Bact – Dithiocarbamates`,
           fung_ben = `Fung & Bact – Benzimidazoles`,
           fung_tri = `Fung & Bact – Triazoles, diazoles`,
           fung_dia = `Fung & Bact – Diazines, morpholines`,
           fung_nes = `Fung & Bact - nes`
           )
  
  
  # Example formulas (replace variable names with yours)
  bf_total <- bf(total_pest ~ ha_cropland_in_2k_radius)                    # Layer 2
  bf_insect <- bf(Insecticides ~ total_pest)                                 # Layer 3 (Insecticides)
  bf_herb   <- bf(Herbicides ~ total_pest)                                   # Layer 3 (Herbicides)
  bf_fung   <- bf(fung_and_bact ~ total_pest)                                   # Layer 3 (Fungicides/Bactericides)
  bf_pgr    <- bf(plant_G ~ total_pest)                                  # Layer 3 (Plant Growth Regulators)
  # # Layer 4 within Insecticides
  # bf_ins_org  <- bf(ins_org ~ Insecticides)
  # bf_ins_car  <- bf(ins_car ~ Insecticides)
  # bf_ins_pyr  <- bf(ins_pyr ~ Insecticides)
  # bf_ins_nes  <- bf(ins_nes ~ Insecticides)
  # # Layer 4 within Herbicides
  # bf_herb_phe  <- bf(herb_phe ~ Herbicides)
  # bf_herb_tri <- bf(herb_tri ~ Herbicides)
  # bf_herb_ami  <- bf(herb_ami ~ Herbicides)
  # bf_herb_car  <- bf(herb_car ~ Herbicides)
  # bf_herb_din  <- bf(herb_din ~ Herbicides)
  # bf_herb_ure  <- bf(herb_ure ~ Herbicides)
  # bf_herb_nes  <- bf(herb_nes ~ Herbicides)
  # # Layer 4 within Fung and Bact
  # bf_fung_ino <- bf(fung_ino ~ fung_and_bact)
  # bf_fung_dit <- bf(fung_dit ~ fung_and_bact)
  # bf_fung_ben <- bf(fung_ben ~ fung_and_bact)
  # bf_fung_tri <- bf(fung_tri ~ fung_and_bact)
  # bf_fung_dia<- bf(fung_dia ~ fung_and_bact)
  # bf_fung_nes <- bf(fung_nes ~ fung_and_bact)
  
  
  # Finally, the outcome model for Ct
  bf_Ct <- bf(Ct | cens(censoring) ~ ha_cropland_in_2k_radius + total_pest +
                Insecticides + Herbicides + fung_and_bact + plant_G + 
                Season + (1 | Hive_ID))
  
  # Combine models into a joint model:
  joint_model[[poi]] <- brm(
    bf_total +
      bf_insect + bf_herb + bf_fung + bf_pgr + 
        # bf_ins_org + bf_ins_car + bf_ins_pyr + bf_ins_nes +
        # bf_herb_phe + bf_herb_tri + bf_herb_ami + bf_herb_car + bf_herb_din + bf_herb_ure + bf_herb_nes +
        # bf_fung_ino + bf_fung_dit + bf_fung_ben + bf_fung_tri + bf_fung_dia + bf_fung_nes +
      bf_Ct + set_rescor(FALSE),
    family = brmsfamily("gaussian", link = "softplus"),
    prior(normal(20, 5), class = "Intercept", lb = 1, ub = 40),
    data = test_tibble_ct[[poi]],        # your data frame containing all these variables
    chains = 4, iter = 10000, cores = 4,
    seed = 1
    )
  
  summary(joint_model[[poi]])
  # plot(joint_model[[poi]])
  # pairs(joint_model[[poi]])
  # pp_check(joint_model[[poi]], resp = "Ct")
  # nuts_params(joint_model[[poi]])
    
  # saveRDS(joint_model[[poi]], paste0("temp.bayes/joint.10k_iter_20frogs.simple.all_layers", poi,".rds"))
  # temp_model <- readRDS(paste0("temp.bayes/joint.20k_iter_10frogs.nestedDWV B.rds"))
  # temp_model <- readRDS(paste0("temp.bayes/joint.10k_iter_10frogs.simpleDWV B.rds"))
  # summary(temp_model)
  # 
    
  # model_ct_censored[[poi]] <- brm(
  #   formula = bf(Ct | cens(censoring) ~ ha_cropland_in_2k_radius + Season + (1 | Hive_ID )),
  #   seed = 1,
  #   data = test_tibble_ct[[poi]],
  #   prior(normal(20, 5), class = "Intercept", lb = 1, ub = 40),
  #   family = brmsfamily("gaussian", link = "softplus"),
  #   # family = gaussian(),
  #   chains = 4,         # 4 independent MCMC chains
  #   iter = 10000,        # 2000 iterations per chain (roughly 1000 warmup, 1000 sampling)
  #   cores = 4           # Utilize 4 CPU cores in parallel
  # )
  

}
