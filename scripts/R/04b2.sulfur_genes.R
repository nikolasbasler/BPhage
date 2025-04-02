library(tidyverse)
library(ggpubr)
library(patchwork)

set.seed(1)

metadata <- readRDS("output/R/R_variables/metadata.RDS")
classification <- readRDS("output/R/R_variables/classification.RDS")
present_in_all_countries <- read_lines("data/core_contigs.txt")

all_hosts <- read.csv("output/R/host_pies/all_hosts.all.csv") %>%
  rename(contig = Virus, Host_genus = Genus)

phold_predictions_with_extensions <- read.csv("output/R/gene_content/phold_predictions_with_extensions.csv")

# phage_tpm <- read.csv("output/R/relative_abundance/phage_tpm.csv")
phage_bee_pool_tpm <- read.csv("output/R/relative_abundance/phage_bee_pool_tpm.csv")

cropland_fraction <- read.csv("data/land_cover_results.csv") %>% 
  tibble() %>%
  mutate(cropland_fraction = cropland_fraction / 100) %>%
  rename(cropland_fraction_2k_radius = cropland_fraction) %>%
  arrange(cropland_fraction_2k_radius)
                            
FAOSTAT_pest_data <- read.csv("data/FAOSTAT_pest_data_en_3-4-2025.csv")
FAOSTAT_area_data <- read.csv("data/FAOSTAT_area_data_en_3-5-2025.csv")
FAOSTAT_pop_data <- read.csv("data/FAOSTAT_pop_data_en_3-5-2025.csv")

landuse <- FAOSTAT_area_data %>% 
  filter(Element %in% c("Area", "Value of agricultural production (Int. $) per Area"), 
         Item %in% c("Land area", "Cropland", "Permanent crops", "Temporary crops", "Agricultural land")) %>%
  filter(!(Element == "Area" & Item == "Agricultural land")) %>%
  mutate(Item = ifelse(Element == "Value of agricultural production (Int. $) per Area", "agri_value", Item)) %>%
  select(Area, Year, Item, Value) %>% 
  pivot_wider(id_cols = c(Area, Year), values_from = Value, names_from = Item) %>%
  mutate(active_crops = `Temporary crops` + `Permanent crops`)
pops <- FAOSTAT_pop_data %>%
  select(Area, Year, Element, Value) %>%
  mutate(Value = Value * 1000) %>%
  pivot_wider(id_cols = c(Area, Year), values_from = Value, names_from = Element)

FAOSTAT_added_data <- FAOSTAT_pest_data %>%
  select(Area, Element, Item, Year, Value) %>%
  filter(!is.na(Value), 
         Year >= 2019) %>%
  pivot_wider(id_cols = c(Area, Item, Year), names_from = Element, values_from = Value) %>%
  left_join(., landuse, by = c("Area", "Year")) %>%
  left_join(., pops, by = c("Area", "Year")) %>%
  mutate(Country = case_when(Area == "Belgium" ~ "BE",
                             Area == "France" ~ "FR",
                             Area == "Germany" ~ "DE",
                             Area == "Netherlands (Kingdom of the)" ~ "NL",
                             Area == "Portugal" ~ "PT",
                             Area == "Romania" ~ "RO",
                             Area == "Switzerland" ~ "CH",
                             Area == "United Kingdom of Great Britain and Northern Ireland" ~ "UK"
                             )) %>%
  mutate(Country = factor(Country, levels = c("PT", "FR", "UK", "BE", "NL", "CH", "DE", "RO"))) %>%
  mutate(`Use per area of cropland` = ifelse(is.na(`Use per area of cropland`), `Agricultural Use` / Cropland, `Use per area of cropland`),
         `Use per capita` = ifelse(is.na(`Use per capita`), `Agricultural Use` / `Total Population - Both sexes`, `Use per capita`),
         `Use per value of agricultural production` = ifelse(is.na(`Use per value of agricultural production`), `Agricultural Use` / agri_value, `Use per value of agricultural production`),
         use_per_active_cropland = `Agricultural Use` / active_crops,
         use_per_land_area = `Agricultural Use` / `Land area`) %>%
  left_join(., cropland_fraction, by = "Country") %>%
  mutate(est_use_in_2k_radius = `Use per area of cropland` * cropland_fraction_2k_radius) %>%
  select(Country, Item, Year, `Agricultural Use`, `Use per area of cropland`, 
         `Use per capita`, `Use per value of agricultural production`, 
         use_per_active_cropland, use_per_land_area, est_use_in_2k_radius)

sulfur_phage_annot <- phold_predictions_with_extensions %>%
  # filter(str_detect(product, "sulf")) %>%
  filter(product == "phosphoadenosine phosphosulfate reductase") %>% 
  inner_join(classification, ., by = join_by("contig" == "contig_id"))

sulf_meta <- phage_bee_pool_tpm %>% tibble() %>%
  pivot_longer(-contig, names_to = "Bee_pool", values_to = "TPM") %>%
  mutate(sulf_metabolism = ifelse(contig %in% sulfur_phage_annot$contig, TRUE, FALSE)) %>%
  left_join(., metadata %>% select(Bee_pool, Hive_ID, Country, Season, Health) %>% distinct(), by = "Bee_pool") %>%
  left_join(., classification, by = "contig") %>%
  select(contig, Bee_pool, TPM, sulf_metabolism, Hive_ID, Country, Season, Health, Family, Host_group, Core) %>%
  # filter(Gut_part == "rec") %>% # Focus on rectum for more stable values (some midgut samples have very few phages)
  group_by(Bee_pool, sulf_metabolism) %>% 
  mutate(sulf_tpm = sum(TPM)) %>%
  mutate(Hive_ID = factor(Hive_ID)) %>%
  ungroup()

meta_variables <- c("Hive_ID", "Country", "Season", "Health")
sulf_stats <-  list()
tpm_KW <- list()
sulf_tpm_plots <- list()
genomes_KW <- list()
sulf_genome_prop_plots <- list()
# nosema_mapped_prop_plots <- list()
for (met_v in meta_variables) {
  tpm <- sulf_meta %>%
    filter(sulf_metabolism) %>%
    group_by(.data[[met_v]]) %>%
    mutate(mean_sulf_tpm = mean(sulf_tpm)) %>%
    filter(mean_sulf_tpm > 0) %>%
    ungroup() %>%
    select(Bee_pool, all_of(meta_variables), sulf_tpm, mean_sulf_tpm) %>%
    distinct() %>%
    arrange(desc(mean_sulf_tpm))
  tpm_KW[[met_v]] <- kruskal.test(tpm$sulf_tpm ~ tpm[[met_v]]) %>%
    unlist() %>% as.matrix() %>% t() %>% as_tibble() %>% 
    rename(chi_squared = `statistic.Kruskal-Wallis chi-squared`) %>%
    mutate(parameter.df = as.numeric(parameter.df),
           p.value = as.numeric(p.value),
           chi_squared = as.numeric(chi_squared))

  genomes <- sulf_meta %>%
    filter(TPM >0) %>%
    group_by(Bee_pool, sulf_metabolism) %>% 
    mutate(sulf_genomes = n()) %>%
    group_by(Bee_pool) %>%
    mutate(sample_sulf_genome_props = sulf_genomes / n ()) %>% 
    group_by(.data[[met_v]]) %>%
    mutate(sulf_genome_count = sum(sulf_metabolism),
           sulf_genome_prop = mean(sulf_metabolism)) %>%
    ungroup() %>%
    filter(sulf_metabolism) %>%
    select(Bee_pool, all_of(meta_variables), sample_sulf_genome_props, sulf_genome_prop) %>%
    distinct() %>%
    arrange(desc(sulf_genome_prop))
  
  genomes_KW[[met_v]] <- kruskal.test(genomes$sample_sulf_genome_props ~ genomes[[met_v]]) %>%
    unlist() %>% as.matrix() %>% t() %>% as_tibble() %>% 
    rename(chi_squared = `statistic.Kruskal-Wallis chi-squared`) %>%
    mutate(parameter.df = as.numeric(parameter.df),
           p.value = as.numeric(p.value),
           chi_squared = as.numeric(chi_squared))

  sulf_tpm_plots[[met_v]] <- tpm %>%
    ggplot(aes(x = .data[[met_v]], y = sulf_tpm)) +
    geom_boxplot(outliers = FALSE) +
    geom_jitter(alpha = 0.33, width = 0.25, height = NULL) +
    ggtitle(paste0("KW p-value: ", tpm_KW[[met_v]]$p.value)) 
  
  sulf_genome_prop_plots[[met_v]] <- genomes %>%
    ggplot(aes(x = .data[[met_v]], y = sample_sulf_genome_props)) +
    geom_boxplot(outliers = FALSE) +
    geom_jitter(alpha = 0.33, width = 0.25) +
    ggtitle(paste0("KW p-value: ", genomes_KW[[met_v]]$p.value))
  
  # Nosema
  # for (set in names(nosema_mapped_prop)) {
  #   nosema_mapped_prop_plots[[met_v]][[set]] <- ggplot(nosema_mapped_prop[[set]], aes(x = .data[[met_v]], y = mapped_to_nosema_prop)) +
  #     geom_boxplot(outliers = FALSE) +
  #     geom_jitter(alpha = 0.33, width = 0.25) +
  #     ggtitle(set)
  # }
  
  sulf_tpm_plots[[paste0(met_v,"_facet")]] <- sulf_tpm_plots[[met_v]] +
    facet_wrap(~Country, scales = "free_x") +
    ggtitle("")
  sulf_genome_prop_plots[[paste0(met_v,"_facet")]] <- sulf_genome_prop_plots[[met_v]]  +
    facet_wrap(~Country, scales = "free_x") +
    ggtitle("")
  # nosema_mapped_prop_plots[[paste0(met_v,"_facet")]][[set]] <- nosema_mapped_prop_plots[[met_v]][[set]] +
  #   facet_wrap(~Country, scales = "free_x") +
  #   ggtitle(set)


  summarised_tpm <- tpm %>%
    select(all_of(met_v), mean_sulf_tpm) %>%
    distinct()
  summarised_genomes <- genomes %>% 
    select(all_of(met_v), sulf_genome_prop) %>% 
    distinct()
  
  sulf_stats[[met_v]] <- left_join(summarised_tpm, summarised_genomes, by = met_v)
 
}

sulf_stats[["Sesaon_facet"]] <- sulf_meta %>%
  filter(sulf_metabolism) %>%
  group_by(Country, Season) %>%
  mutate(mean_sulf_season_tpm = mean(sulf_tpm)) %>%
  filter(mean_sulf_season_tpm > 0) %>%
  ungroup() %>%
  select(Country, Season, mean_sulf_season_tpm) %>%
  distinct()
total_genomes <- sulf_meta %>%
  filter(TPM > 0) %>%
  group_by(Country, Season) %>%
  summarise(total_genome_count = n(), .groups = "drop")
sulf_stats[["Sesaon_facet"]] <- sulf_meta %>%
  filter(TPM > 0,
         sulf_metabolism) %>%
  group_by(Country, Season) %>%
  summarise(genome_count = n(), .groups = "drop") %>%
  full_join(., total_genomes, by = c("Country", "Season")) %>%
  mutate(sulf_genome_prop = genome_count / total_genome_count) %>%
  select(Country, Season, sulf_genome_prop) %>%
  left_join(sulf_stats[["Sesaon_facet"]], ., by = c("Country", "Season"))
  

sulf_stats_plots <- list()
sulf_stats_cor_plots <- list()
for(var in c("Hive_ID", "Country", "Season", "Health")) {
  sulf_stats_plots[[var]] <- sulf_stats[[var]] %>%
    pivot_longer(-all_of(var), names_to = "measure", values_to = "mean_tpm") %>%
    ggplot(aes(x = .data[[var]], y = mean_tpm, fill = measure)) +
    geom_col(position = "dodge")
  
  sulf_stats_cor_plots[[var]] <- sulf_stats[[var]] %>%
    ggplot(aes(x = mean_sulf_tpm,  y = sulf_genome_prop)) +
    geom_point() +
    geom_smooth(method = "glm", formula = y ~ x) +
    stat_cor(method = "spearman") +
    ggtitle(var)

}

sulf_stats_plots[["Season_facet"]] <- sulf_stats$Sesaon_facet %>%
  pivot_longer(-c(Country, Season), names_to = "measure", values_to = "mean_tpm") %>%
  ggplot(aes(x = Season, y = mean_tpm, fill = measure)) +
  geom_col(position = "dodge") +
  facet_wrap(~Country, scales = "free_x")

sulf_stats$Country %>%
  full_join(., cropland_fraction, by = "Country") %>%
  ggplot(aes(x = cropland_fraction_2k_radius, y = mean_sulf_tpm)) +
  geom_point() +
  geom_smooth(method = "glm", formula = y ~ x) +
  stat_cor(method = "spearman")


## Pesticide data

elements <- c("Agricultural Use", "Use per area of cropland", "Use per capita",
              "Use per value of agricultural production", 
              "use_per_active_cropland", "use_per_land_area", "est_use_in_2k_radius")
# elements <- c("Use per area of cropland", "Use per capita")
all_insecticides <- FAOSTAT_added_data %>% 
  distinct(Item) %>% 
  filter(str_detect(Item, "Insecticides")) %>% 
  unlist(use.names = FALSE)

FAOSTAT_plots <- list()
for (element in elements) {
  FAOSTAT_plots[[element]] <- FAOSTAT_added_data %>% 
    select(all_of(c("Country", "Year", "Item", element))) %>%
    # filter(Year >= 2019,
    #        Item %in% all_insecticides) %>%
    filter(Year == 2019) %>%
    mutate(Year = as.character(Year)) %>%
    ggplot(aes(x = Country, y = .data[[element]], fill = Year)) +
    geom_col(position = "dodge") +
    ggtitle(element) +
    facet_wrap(~Item, scales = "free")
}



pest_cor_raw_p <- list()
pest_cor_tibbles <-list()
pest_cor_plots <- list()
pest_cor_tests <- list()
elements <- c("Use per capita", "Use per area of cropland", "use_per_land_area", "est_use_in_2k_radius")
for (element in elements) {
  # for (item in unique(FAOSTAT_added_data$Item)) {
  # for (item in c("Insecticides", "Herbicides", "Fungicides and Bactericides")) {
  for (item in c("Pesticides (total)", "Insecticides", "Herbicides", "Fungicides and Bactericides")) {
    # for (item in c("Insecticides")) {
    # for (item in all_insecticides) {
    # for (year in as.character(2017:2020)) {
    # for (year in as.character(2019:2020)) {
    for (year in as.character(2019)) {
      cor_tibble <- FAOSTAT_added_data %>%
        select(all_of(c("Country", "Year", "Item", element))) %>%
        filter(Item == item) %>%
        filter(Year == year) %>%
        left_join(., sulf_stats$Country, by = "Country") %>%
        select(-sulf_genome_prop)
      
      if (nrow(cor_tibble) > 1) {
        pest_cor_raw_p[[paste0(element, "/", item, "/", year)]] <- cor.test(cor_tibble[[element]], cor_tibble$mean_sulf_tpm, method = "spearman")$p.value
        # pest_cor_raw_p[[element]][[item]][[year]] <- cor.test(cor_tibble[[element]], cor_tibble$mean_sulf_tpm, method = "spearman")$p.value
        correlation <- cor.test(cor_tibble[[element]], cor_tibble$mean_sulf_tpm, method = "spearman")
        # if (correlation$p.value <= 0.05) {
        pest_cor_tests[[element]][[item]][[year]] <- correlation
        pest_cor_tibbles[[element]][[item]][[year]] <- cor_tibble
        pest_cor_plots[[element]][[item]][[year]] <- cor_tibble %>%
          ggplot(aes(x = .data[[element]], y = mean_sulf_tpm)) +
          geom_point() +
          geom_smooth(method = "glm", formula = y ~ x) +
          stat_cor(method = "spearman", label.x = 0, label.y = max(cor_tibble$mean_sulf_tpm)*1.1) +
          ggtitle(paste0(element,  " - ", item, " - ", year))
      }
    }
  }
}

narrowed_data <- FAOSTAT_added_data %>% 
  # filter(Year == 2020) %>%
  filter(Year == 2019) %>%
  # filter(Item %in% c("Pesticides (total)", "Insecticides", "Herbicides", "Fungicides and Bactericides"))
  # filter(Item %in% c("Pesticides (total)", "Herbicides", "Fungicides and Bactericides") | str_detect(Item, "Insecticides"))
  filter(Item %in% unique(FAOSTAT_added_data$Item))

pest_cor_season_tibble <- sulf_meta %>%
  filter(sulf_metabolism) %>%
  group_by(Country, Season) %>%
  mutate(c_s_mean_sulf_tpm = mean(sulf_tpm)) %>%
  filter(c_s_mean_sulf_tpm > 0) %>%
  ungroup() %>%
  select(Country, Season, c_s_mean_sulf_tpm) %>%
  distinct() %>%
  arrange(desc(c_s_mean_sulf_tpm)) %>%
  full_join(., narrowed_data, by = "Country", relationship = "many-to-many") 


season_cor_plots <- list()
season_cor_tests <- list()
for (item in unique(pest_cor_season_tibble$Item)) {
  # for (element in c("Use per area of cropland", "Use per capita", "Use per value of agricultural production", "use_per_active_cropland", "use_per_land_area")) {
  # for (element in c("Use per capita", "Use per area of cropland", "use_per_active_cropland", "use_per_land_area")) {
  # for (element in c("Use per capita", "Use per area of cropland", "use_per_land_area")) {
  # for (element in c("use_per_land_area", "est_use_in_2k_radius")) {
  for (element in c("est_use_in_2k_radius")) {
    plot_tibble <- pest_cor_season_tibble %>%
      filter(Item == item) %>%
      select(Country, Season, c_s_mean_sulf_tpm, all_of(element))
    
    season_cor_plots[[item]][[element]] <- plot_tibble %>%
      ggplot(aes(x = .data[[element]], y = c_s_mean_sulf_tpm)) +
      geom_point() +
      facet_wrap(~Season, scales = "free") +
      geom_smooth(method = "glm", formula = y ~ x) +
      stat_cor(method = "spearman") +
      ggtitle(paste0(item, " - ", element))
    
    # for (season in c("spr", "sum", "aut")) {
    #   cor_tibble <- plot_tibble %>%
    #     filter(Season == season)
    #   if (nrow(cor_tibble) > 1) {
    #     season_cor_tests[[item]][[element]][[season]] <- cor.test(cor_tibble$c_s_mean_sulf_tpm, cor_tibble$use_per_land_area, method = "spearman")$p.value
    #   }
    # }
  }
}

