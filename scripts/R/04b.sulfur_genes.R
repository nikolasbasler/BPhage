library(tidyverse)
library(ggpubr)
library(multcompView)
library(patchwork)

set.seed(1)

metadata <- readRDS("output/R/R_variables/metadata.RDS")
classification <- readRDS("output/R/R_variables/classification.RDS")
present_in_all_countries <- read_lines("data/core_contigs.txt")

all_hosts <- read.csv("output/R/host_pies/all_hosts.all.csv") %>%
  rename(contig = Virus, Host_genus = Genus)

phold_predictions_with_extensions <- read.csv("output/R/gene_content/phold_predictions_with_extensions.csv")

phage_tpm <- read.csv("output/R/relative_abundance/phage_tpm.csv")

presence_absence <- list()
presence_absence$Countries <- read.csv("output/R/prevalence/presence_absence.Countries.csv")
presence_absence$Seasons <- read.csv("output/R/prevalence/presence_absence.Seasons.csv")

cropland_fraction <- read.csv("data/land_cover_results.csv") %>% 
  tibble() %>%
  mutate(cropland_fraction = cropland_fraction / 100) %>%
  rename(cropland_fraction_2k_radius = cropland_fraction)

FAOSTAT_pest_data <- read.csv("data/FAOSTAT_pest_data_en_3-4-2025.csv") %>%
  mutate(Country = case_when(Area == "Belgium" ~ "BE",
                          Area == "France" ~ "FR",
                          Area == "Germany" ~ "DE",
                          Area == "Netherlands (Kingdom of the)" ~ "NL",
                          Area == "Portugal" ~ "PT",
                          Area == "Romania" ~ "RO",
                          Area == "Switzerland" ~ "CH",
                          Area == "United Kingdom of Great Britain and Northern Ireland" ~ "UK"
                          )
         ) %>%
  mutate(Country = factor(Country, levels = c("PT", "FR", "UK", "BE", "NL", "CH", "DE", "RO"))) %>%
  tibble()
FAOSTAT_area_data <- read.csv("data/FAOSTAT_area_data_en_3-5-2025.csv") %>%
  mutate(Country = case_when(Area == "Belgium" ~ "BE",
                             Area == "France" ~ "FR",
                             Area == "Germany" ~ "DE",
                             Area == "Netherlands (Kingdom of the)" ~ "NL",
                             Area == "Portugal" ~ "PT",
                             Area == "Romania" ~ "RO",
                             Area == "Switzerland" ~ "CH",
                             Area == "United Kingdom of Great Britain and Northern Ireland" ~ "UK"
  )
  ) %>%
  mutate(Country = factor(Country, levels = c("PT", "FR", "UK", "BE", "NL", "CH", "DE", "RO"))) %>%
  tibble()
FAOSTAT_pop_data <- read.csv("data/FAOSTAT_pop_data_en_3-5-2025.csv") %>%
  mutate(Country = case_when(Area == "Belgium" ~ "BE",
                             Area == "France" ~ "FR",
                             Area == "Germany" ~ "DE",
                             Area == "Netherlands (Kingdom of the)" ~ "NL",
                             Area == "Portugal" ~ "PT",
                             Area == "Romania" ~ "RO",
                             Area == "Switzerland" ~ "CH",
                             Area == "United Kingdom of Great Britain and Northern Ireland" ~ "UK"
  )
  ) %>%
  mutate(Country = factor(Country, levels = c("PT", "FR", "UK", "BE", "NL", "CH", "DE", "RO"))) %>%
  tibble()
landuse <- FAOSTAT_area_data %>% 
  filter(Element %in% c("Area", "Value of agricultural production (Int. $) per Area"), 
         Item %in% c("Land area", "Cropland", "Permanent crops", "Temporary crops", "Agricultural land")) %>%
  filter(!(Element == "Area" & Item == "Agricultural land")) %>%
  mutate(Item = ifelse(Element == "Value of agricultural production (Int. $) per Area", "agri_value", Item)) %>%
  select(Country, Year, Item, Value) %>% 
  pivot_wider(id_cols = c(Country, Year), values_from = Value, names_from = Item) %>%
  mutate(active_crops = `Temporary crops` + `Permanent crops`)
pops <- FAOSTAT_pop_data %>%
  select(Country, Year, Element, Value) %>%
  mutate(Value = Value * 1000) %>%
  pivot_wider(id_cols = c(Country, Year), values_from = Value, names_from = Element)

# cropland <- FAOSTAT_pest_data %>% 
#   filter(Element %in% c("Agricultural Use", "Use per area of cropland"),
#          Item == "Pesticides (total)") %>%
#   pivot_wider(id_cols = c(Country, Year), names_from = Element, values_from = Value) %>%
#   mutate(cropland = `Agricultural Use` * 1000 / `Use per area of cropland`) %>%
#   select(Country, Year, cropland)
# populations <- FAOSTAT_pest_data %>% 
#   filter(Element %in% c("Agricultural Use", "Use per capita"),
#          Item == "Pesticides (total)") %>%
#   pivot_wider(id_cols = c(Country, Year), names_from = Element, values_from = Value) %>%
#   mutate(population = `Agricultural Use` * 1000 / `Use per capita`) %>%
#   select(Country, Year, population)
# agri_values <- FAOSTAT_pest_data %>% 
#   filter(Element %in% c("Agricultural Use", "Use per value of agricultural production"),
#          Item == "Pesticides (total)") %>%
#   pivot_wider(id_cols = c(Country, Year), names_from = Element, values_from = Value) %>%
#   mutate(agri_values = `Agricultural Use` / `Use per value of agricultural production`) %>%
#   select(Country, Year, agri_values)

# FAOSTAT_added_data <- FAOSTAT_pest_data %>% 
#   select(Country, Element, Item, Year, Value) %>%
#   filter(!is.na(Value)) %>%
#   pivot_wider(id_cols = c(Country, Item, Year), names_from = Element, values_from = Value) %>%
#   left_join(., populations, by = c("Country", "Year")) %>%
#   left_join(., cropland, by = c("Country", "Year")) %>%
#   left_join(., agri_values, by = c("Country", "Year")) %>%
#   mutate(`Use per area of cropland` = ifelse(is.na(`Use per area of cropland`), `Agricultural Use` / cropland, `Use per area of cropland`),
#          `Use per capita` = ifelse(is.na(`Use per capita`), `Agricultural Use` / population, `Use per capita`),
#          `Use per value of agricultural production` = ifelse(is.na(`Use per value of agricultural production`), `Agricultural Use` / agri_values, `Use per value of agricultural production`)
#          )

FAOSTAT_added_data <- FAOSTAT_pest_data %>%
  select(Country, Element, Item, Year, Value) %>%
  filter(!is.na(Value), 
         Year >= 2019) %>%
  pivot_wider(id_cols = c(Country, Item, Year), names_from = Element, values_from = Value) %>%
  left_join(., landuse, by = c("Country", "Year")) %>%
  left_join(., pops, by = c("Country", "Year")) %>%
  mutate(`Use per area of cropland` = ifelse(is.na(`Use per area of cropland`), `Agricultural Use` / Cropland, `Use per area of cropland`),
         `Use per capita` = ifelse(is.na(`Use per capita`), `Agricultural Use` / `Total Population - Both sexes`, `Use per capita`),
         `Use per value of agricultural production` = ifelse(is.na(`Use per value of agricultural production`), `Agricultural Use` / agri_value, `Use per value of agricultural production`),
         use_per_active_cropland = `Agricultural Use` / active_crops,
         use_per_land_area = `Agricultural Use` / `Land area`
         ) %>%
  select(Country, Item, Year, `Agricultural Use`, `Use per area of cropland`, 
         `Use per capita`, `Use per value of agricultural production`, use_per_active_cropland, use_per_land_area)

# geom_col()# Nosema
nosema_mapped_counts <- list()
nosema_mapped_counts$pools <- read.delim("output/nosema_mapped_counts.tsv") %>% tibble()
nosema_mapped_counts$rec <- read.delim("output/nosema_mapped_counts_rec.tsv") %>% tibble()

nosema_mapped_prop <- list()
nosema_mapped_prop$pools <- left_join(nosema_mapped_counts$pools, metadata, by = "Bee_pool") %>%
  mutate(mapped_to_nosema_prop = mapped_to_Nosema / Hostout_R1_plus_R2) %>%
  mutate(Hive_ID = factor(Hive_ID)) %>%
  select(Bee_pool, mapped_to_nosema_prop, Hive_ID, Country, Season, Health) %>%
  distinct()
nosema_mapped_prop$rec <- left_join(nosema_mapped_counts$rec, metadata, by = "Sample_ID") %>%
  mutate(mapped_to_nosema_prop = mapped_to_Nosema / Hostout_R1_plus_R2) %>%
  mutate(Hive_ID = factor(Hive_ID)) %>%
  select(Sample_ID, Bee_pool, mapped_to_nosema_prop, Hive_ID, Country, Season, Health) %>%
  distinct()

sulfur_phage_annot <- phold_predictions_with_extensions %>%
  # filter(str_detect(product, "sulf")) %>%
  filter(product == "phosphoadenosine phosphosulfate reductase") %>% 
  inner_join(classification, ., by = join_by("contig" == "contig_id"))

sulf_meta <- phage_tpm %>%
  pivot_longer(-contig, names_to = "Sample_ID", values_to = "TPM") %>%
  mutate(sulf_metabolism = ifelse(contig %in% sulfur_phage_annot$contig, TRUE, FALSE)) %>%
  left_join(., metadata, by = "Sample_ID") %>% 
  select(contig, Sample_ID, Bee_pool, TPM, sulf_metabolism, Hive_ID, Country, Season, Gut_part, Health) %>%
  filter(Gut_part == "rec") %>% # Focus on rectum for more stable values (some midgut samples have very few phages)
  group_by(Sample_ID, sulf_metabolism) %>% 
  mutate(sulf_tpm = sum(TPM)) %>%
  mutate(Hive_ID = factor(Hive_ID)) %>%
  ungroup()

nosema_sulf_cor_tibble <- list()
nosema_sulf_cor_tibble$pools <- sulf_meta %>% 
  filter(sulf_metabolism) %>%
  group_by(Bee_pool) %>%
  mutate(sulf_bee_pool_tpm = mean(sulf_tpm)) %>%
  ungroup() %>%
  select(Bee_pool, sulf_bee_pool_tpm) %>%
  distinct() %>%
  left_join(., nosema_mapped_prop$pools, by = "Bee_pool") %>%
  select(Bee_pool, sulf_bee_pool_tpm, mapped_to_nosema_prop) %>%
  distinct()
  
nosema_sulf_cor_tibble$rec <- full_join(nosema_mapped_prop$rec, sulf_meta, by = "Bee_pool") %>%
  filter(Gut_part == "rec") %>% # Focus on rectum for more stable values (some midgut samples have very few phages)
  filter(sulf_metabolism) %>%
  select(Bee_pool, sulf_tpm, mapped_to_nosema_prop) %>%
  distinct()
# nosema_sulf_cor_tibble <- full_join(nosema_mapped_prop, sulf_meta, by = "Sample_ID") %>%
#   filter(sulf_metabolism) %>%
#   select(Sample_ID, sulf_tpm, mapped_to_nosema_prop) %>%
#   distinct()

nosema_sulf_cor_plot <- list()
nosema_sulf_cor_plot$rec <- nosema_sulf_cor_tibble$rec %>%
  ggplot(aes(x = mapped_to_nosema_prop, y = sulf_tpm)) +
  geom_point() +
  stat_cor(method = "spearman", label.x = 0, label.y = max(nosema_sulf_cor_tibble$rec$sulf_tpm)) +  # Add correlation coefficient
  geom_smooth(method = "glm", formula = y ~ x) +
  labs(x = "Proportion of reads mapping to Varimorpha genome", y = "Mean relative abundance of sulfur phages", title = "rec")
nosema_sulf_cor_plot$pools <- nosema_sulf_cor_tibble$pools %>%
  ggplot(aes(x = mapped_to_nosema_prop, y = sulf_bee_pool_tpm)) +
  geom_point() +
  stat_cor(method = "spearman", label.x = 0, label.y = max(nosema_sulf_cor_tibble$pools$sulf_bee_pool_tpm)) +  # Add correlation coefficient
  geom_smooth(method = "glm", formula = y ~ x) +
  labs(x = "Proportion of reads mapping to Varimorpha genome", y = "Mean relative abundance of sulfur phages", title = "pools")

sulf_positive_hives <- sulf_meta %>%
  filter(sulf_tpm > 0) %>%
  group_by(Country) %>%
  summarise(hive_count = n_distinct(Hive_ID),
            sulf_positive_hives = n_distinct(Hive_ID[sulf_metabolism])) %>%
  mutate(sulf_positive_prop = sulf_positive_hives / hive_count)

meta_variables <- c("Hive_ID", "Country", "Season", "Health")
sulf_stats <-  list()
tpm_KW <- list()
tpm_pwc <- list()
sulf_tpm_plots <- list()
genomes_KW <- list()
genomes_pwc <- list()
sulf_genome_prop_plots <- list()
nosema_mapped_prop_plots <- list()
for (met_v in meta_variables) {
  tpm <- sulf_meta %>%
    filter(sulf_metabolism) %>%
    group_by(.data[[met_v]]) %>%
    mutate(mean_sulf_tpm = mean(sulf_tpm)) %>%
    filter(mean_sulf_tpm > 0) %>%
    ungroup() %>%
    select(Sample_ID, all_of(meta_variables), sulf_tpm, mean_sulf_tpm) %>%
    distinct() %>%
    arrange(desc(mean_sulf_tpm))
  tpm_KW[[met_v]] <- kruskal.test(tpm$sulf_tpm ~ tpm[[met_v]]) %>%
    unlist() %>% as.matrix() %>% t() %>% as_tibble() %>% 
    rename(chi_squared = `statistic.Kruskal-Wallis chi-squared`) %>%
    mutate(parameter.df = as.numeric(parameter.df),
           p.value = as.numeric(p.value),
           chi_squared = as.numeric(chi_squared))
  
  # Unfortunately, this still doesn't work. I don't understand how the input for multcompLetters()
  # needs to be. When feeding it the output of pairwise.wilcox.test()$p.value it ignores one of 
  # the categories When feeding it a complete symmetric df, it puts all categories into group "a"...
  # Maybe this helps: https://www.r-bloggers.com/2017/03/perform-pairwise-wilcoxon-test-classify-groups-by-significance-and-plot-results/
  # pwc_p <- pairwise.wilcox.test(tpm$sulf_tpm, tpm[[met_v]], p.adjust.method = "BH")$p.value
  # multcompLetters(pwc_p)
  # complete_triangle <- function(matr) {
  #   missing_col <- rownames(matr) %>% 
  #     tail(n=1)
  #   missing_row <- colnames(matr) %>% 
  #     head(n=1)
  #   full <- matr %>%
  #     as.data.frame() %>%
  #     mutate(!!missing_col := NA) %>%
  #     t() %>%
  #     as.data.frame() %>%
  #     mutate(!!missing_row := NA) %>%
  #     relocate(all_of(missing_row), .before = everything()) %>%
  #     t()
  #   full[upper.tri(full)] <- t(full)[upper.tri(full)]
  #   return(full)
  # }
  # pwc_full <- complete_triangle(pwc_p)
  # genomes_pwc[[met_v]] <- multcompLetters(pwc_p_full)$Letters
  
  genomes <- sulf_meta %>%
    filter(TPM >0) %>%
    group_by(Sample_ID, sulf_metabolism) %>% 
    mutate(sulf_genomes = n()) %>%
    group_by(Sample_ID) %>%
    mutate(sample_sulf_genome_props = sulf_genomes / n ()) %>% 
    group_by(.data[[met_v]]) %>%
    mutate(sulf_genome_count = sum(sulf_metabolism),
           sulf_genome_prop = mean(sulf_metabolism)) %>%
    ungroup() %>%
    filter(sulf_metabolism) %>%
    select(Sample_ID, all_of(meta_variables), sample_sulf_genome_props, sulf_genome_prop) %>%
    distinct() %>%
    arrange(desc(sulf_genome_prop))
  genomes_KW[[met_v]] <- kruskal.test(genomes$sample_sulf_genome_props ~ genomes[[met_v]]) %>%
    unlist() %>% as.matrix() %>% t() %>% as_tibble() %>% 
    rename(chi_squared = `statistic.Kruskal-Wallis chi-squared`) %>%
    mutate(parameter.df = as.numeric(parameter.df),
           p.value = as.numeric(p.value),
           chi_squared = as.numeric(chi_squared))
  # pwc_p <- pairwise.wilcox.test(genomes$sample_sulf_genome_props, genomes[[met_v]], p.adjust.method = "BH")$p.value
  # pwc_full <- complete_triangle(pwc_p)
  # genomes_pwc[[met_v]] <- multcompLetters(pwc_p_full)$Letters
  
  sulf_tpm_plots[[met_v]] <- tpm %>%
    # mutate(!!met_v := factor(.data[[met_v]], levels = unique(tpm[[met_v]]))) %>%
    ggplot(aes(x = .data[[met_v]], y = sulf_tpm)) +
    geom_boxplot(outliers = FALSE) +
    geom_jitter(alpha = 0.33, width = 0.25) +
    # scale_y_continuous(trans='log10') + # This would change the picture because 0s would be ignored.
    ggtitle(paste0("KW p-value: ", tpm_KW[[met_v]]$p.value)) # +
  # scale_y_continuous(limits = c(0,0.1))
  
  sulf_genome_prop_plots[[met_v]] <- genomes %>%
    # mutate(!!met_v := factor(.data[[met_v]], levels = unique(genomes[[met_v]]))) %>%
    ggplot(aes(x = .data[[met_v]], y = sample_sulf_genome_props)) +
    geom_boxplot(outliers = FALSE) +
    geom_jitter(alpha = 0.33, width = 0.25) +
    ggtitle(paste0("KW p-value: ", genomes_KW[[met_v]]$p.value))
  
  # Nosema
  for (set in names(nosema_mapped_prop)) {
    nosema_mapped_prop_plots[[met_v]][[set]] <- ggplot(nosema_mapped_prop[[set]], aes(x = .data[[met_v]], y = mapped_to_nosema_prop)) +
      geom_boxplot(outliers = FALSE) +
      geom_jitter(alpha = 0.33, width = 0.25) +
      ggtitle(set)
  }
  
  if (met_v %in% c("Hive_ID", "Season")) {
    sulf_tpm_plots[[paste0(met_v,"_facet")]] <- sulf_tpm_plots[[met_v]] +
      facet_wrap(~Country, scales = "free_x") +
      ggtitle("")
    sulf_genome_prop_plots[[paste0(met_v,"_facet")]] <- sulf_genome_prop_plots[[met_v]]  +
      facet_wrap(~Country, scales = "free_x") +
      ggtitle("")
    for (set in names(nosema_mapped_prop)) {
      nosema_mapped_prop_plots[[paste0(met_v,"_facet")]][[set]] <- nosema_mapped_prop_plots[[met_v]][[set]] +
        facet_wrap(~Country, scales = "free_x") +
        ggtitle(set)
    }
  }
  
  summarised_tpm <- tpm %>%
    select(all_of(met_v), mean_sulf_tpm) %>%
    distinct()
  summarised_genomes <- genomes %>% 
    select(all_of(met_v), sulf_genome_prop) %>% 
    distinct()
  sulf_stats[[met_v]] <- left_join(summarised_tpm, summarised_genomes, by = met_v)
}


sulf_meta_and_host <- classification %>% tibble() %>%
  left_join(., all_hosts, by = "contig") %>%
  mutate(sulf_metabolism = ifelse(contig %in% sulfur_phage_annot$contig, TRUE, FALSE)) %>%
  select(contig, Host_genus, sulf_metabolism)

all_genomes_per_host <- sulf_meta_and_host %>%
  group_by(Host_genus) %>%
  summarise(total_genome_count = n()) %>%
  ungroup() %>%
  mutate(total_genome_prop = total_genome_count / sum(total_genome_count))

sulf_genomes_per_host <- sulf_meta_and_host %>%
  filter(sulf_metabolism) %>%
  group_by(Host_genus) %>%
  summarise(sulf_genome_count = n()) %>%
  mutate(sulf_genome_prop = sulf_genome_count / sum(sulf_genome_count))


full_join(all_genomes_per_host, sulf_genomes_per_host, by = "Host_genus") %>%
  mutate(sulf_genome_count = ifelse(is.na(sulf_genome_count), 0, sulf_genome_count),
         sulf_genome_prop = ifelse(is.na(sulf_genome_prop), 0, sulf_genome_prop)) %>%
  # select(Host_genus, total_genome_prop, sulf_genome_prop) %>%
  mutate(ratio = sulf_genome_prop / total_genome_prop) 
  
# 

contigency_table <- classification %>% tibble() %>%
  mutate(sulf_metabolism = ifelse(contig %in% sulfur_phage_annot$contig, "yes", "no")) %>%
  left_join(., all_hosts, by ="contig") %>% 
  select(contig, Host_genus, sulf_metabolism) %>%
  filter(Host_genus != "unknown") %>%
  group_by(Host_genus, sulf_metabolism) %>%
  summarise(counts = n(), .groups = "drop") %>%
  # mutate(sulf_metabolism = ifelse(sulf_metabolism, "yes", "no")) %>%
  pivot_wider(names_from = sulf_metabolism, values_from = counts, values_fill = 0) %>% 
  column_to_rownames("Host_genus") %>% 
  as.matrix()

ftest <- fisher.test(contigency_table, simulate.p.value = TRUE)
prop.table(contigency_table, margin = 2)

cont_filt <- contigency_table %>%
  as.data.frame() %>%
  filter(yes >4, no >4) %>%
  as.matrix()

csq <- chisq.test(cont_filt)
csq$stdres
csq$residuals

# for_association_test <- phold_predictions_with_extensions %>%
#   tibble() %>%
#   filter(str_starts(contig_id, "NODE")) %>%
#   rename(contig = contig_id) %>%
#   left_join(., classification, by ="contig") %>%
#   select(contig, product, Family, Core) %>%
#   left_join(., all_hosts, by ="contig") %>%
#   left_join(., presence_absence$Countries, by = "contig") %>%
#   left_join(., presence_absence$Seasons, by = "contig") %>%
#   mutate(sulf_metabolism = ifelse(contig %in% sulfur_phage_annot$contig, TRUE, FALSE))
#   

meta_variables <- c("Country", "Season", "Host_genus", "Core")

for (met_v in meta_variables) {
  for_association_test %>%
    select(contig, all_of(met_v), sulf_metabolism) %>%
    # distinct() %>%
    filter(.data[[met_v]] != "unknown") %>%
    group_by(.data[[met_v]], sulf_metabolism) %>%
    summarise(count = sum(sulf_metabolism)) %>% View()
    
}


# for (goi in unique(for_glm$product)) {
# for (goi in "phosphoadenosine phosphosulfate reductase") {
  # for_association_test_filt <- for_association_test %>%
  #   select(contig, Host_genus) %>%
  #   filter(Host_genus != "unknown") %>%
  #   mutate(Host_genus = as.factor(Host_genus)) %>%
  #   mutate(sulfur_metabolism = ifelse(contig %in% sulfur_phage_annot$contig, TRUE, FALSE)) %>%
  #   distinct()
  # 
  # model <- glm(sulfur_metabolism ~ Host_genus, data = for_glm_filt, family = binomial(link = "logit"))
  # summary(model)
  # exp(coef(model))
# }

### Toxins. Not much here, I think.
# 
# 
# phold_predictions_with_extensions %>% tibble() %>%
#   filter(str_starts(contig_id, "NODE")) %>%
#   filter(str_detect(product, "toxin")) %>%
#   select(contig_id, product) %>%
#   distinct() %>%
#   rename(contig = contig_id) %>%
#   left_join(., all_hosts, by = "contig") %>%
#   group_by(product, Genus) %>%
#   summarise(genome_count = n(), .groups = "drop") %>%
#   arrange(desc(product)) -> hosts_of_toxin_phages

## Pesticide data

elements <- c("Agricultural Use", "Use per area of cropland", "Use per capita", "Use per value of agricultural production", "use_per_active_cropland", "use_per_land_area")
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

# LONG STORY SHORT: CORRELATIONS PRETTY MUCH ONLY WITH "Insecticides" (I.E. NO SPECIFICT ONE) FROM 2019 and 2020
# Chlorinated Hydrocarbons correlate strongly but only in "Agricultural Use", i.e. not corrected for population
# or crop land. Also, for 2019 two countries have no data on the usage, for 2020 one country. So I didn't follow it 
# further.

pest_cor_raw_p <- list()
pest_cor_tibbles <-list()
pest_cor_plots <- list()
pest_cor_tests <- list()
# elements <- c("Use per area of cropland", "Use per capita")
# elements <- c("Agricultural Use", "Use per capita", "Use per area of cropland", "Use per value of agricultural production", "use_per_active_cropland", "use_per_land_area")
elements <- c("Use per capita", "Use per area of cropland", "use_per_land_area")
# elements <- c("use_per_land_area")
for (element in elements) {
  # for (item in unique(FAOSTAT_added_data$Item)) {
  # for (item in c("Insecticides", "Herbicides", "Fungicides and Bactericides")) {
  for (item in c("Insecticides", "Herbicides", "Fungicides and Bactericides")) {
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
        # }
      }
    }
  }
}

tibble(element = names(pest_cor_raw_p), p_value = unlist(pest_cor_raw_p)) %>%
  mutate(p_adjust = p.adjust(p_value, method = "BH")) %>% 
  arrange(p_value) %>% 
  filter(p_adjust <= 0.05)

pest_cor_tibble <- FAOSTAT_added_data %>%
  filter(Year >= 2019,
         Item == "Insecticides") %>%
  select(all_of(c("Country", "Year", "Use per capita"))) %>%
  left_join(., sulf_stats$Country, by = "Country") %>%
  select(-sulf_genome_prop) %>%
  pivot_wider(id_cols = c("Country", "mean_sulf_tpm"), names_from = Year, values_from = `Use per capita`) %>%
  rename(use_per_cap_2019 = `2019`,
         use_per_cap_2020 = `2020`)
pest_cor_test_2019 <- cor.test(pest_cor_tibble$mean_sulf_tpm, pest_cor_tibble$use_per_cap_2019)
pest_cor_test_2020 <- cor.test(pest_cor_tibble$mean_sulf_tpm, pest_cor_tibble$use_per_cap_2020)

pest_cor_plot <- pest_cor_tibble %>%
  ggplot() +
  geom_point(aes(x = use_per_cap_2019, y = mean_sulf_tpm, color = "2019")) +
  geom_smooth(aes(x = use_per_cap_2019, y = mean_sulf_tpm, color = "2019"), 
              method = "glm", formula = y ~ x) +
  geom_point(aes(x = use_per_cap_2020, y = mean_sulf_tpm, color = "2020")) +
  geom_smooth(aes(x = use_per_cap_2020, y = mean_sulf_tpm, color = "2020"), 
              method = "glm", formula = y ~ x) +
  stat_cor(aes(x = use_per_cap_2019, y = mean_sulf_tpm), method = "spearman", 
           label.x = 0, label.y = max(pest_cor_tibble$mean_sulf_tpm) * 1.1, color = "blue") +
  stat_cor(aes(x = use_per_cap_2020, y = mean_sulf_tpm), method = "spearman", 
           label.x = 0, label.y = max(pest_cor_tibble$mean_sulf_tpm), color = "red") +
  scale_color_manual(values = c("2019" = "blue", "2020" = "red"), name = "Year") +
  labs(x = "Insecticide use per capita",
       y = "Mean relative abundance of sulfur phages") +
  theme_minimal()

narrowed_data <- FAOSTAT_added_data %>% 
  filter(Year == 2020) %>%
  # filter(Year == 2019) %>%
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
  for (element in c("use_per_land_area")) {
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
      
      for (season in c("spr", "sum", "aut")) {
        cor_tibble <- plot_tibble %>%
          filter(Season == season)
        if (nrow(cor_tibble) > 1) {
          season_cor_tests[[item]][[element]][[season]] <- cor.test(cor_tibble$c_s_mean_sulf_tpm, cor_tibble$use_per_land_area, method = "spearman")$p.value
        }
      }
  }
}

wraps <- list()
for (item in names(season_cor_plots)) {
  wraps[[item]] <- wrap_plots(season_cor_plots[[item]], ncol = 1)
}

# Fung & Bact – Benzimidazoles in summer!
# Herbicides – Triazines in autumn
# Insecticides - nes in autum
# Insecticides - Biopesticides remarkable flat line but mainly driven by 1 country
# Insecticides – Organo-phosphates in summer but mainly driven by 1 country
# Insecticdes (all) strongly correlated in autumn

# Procedure
# 1. All pesticides together
wraps$`Pesticides (total)`



tibble(season = c("spr", "sum", "aut"),
       p_value =  season_cor_tests$`Pesticides (total)` %>% unlist(),
       p_adj = p.adjust(season_cor_tests$`Pesticides (total)` %>% unlist(), method = "BH"))

# 2. Separate groups
# wraps$Insecticides
# wraps$Herbicides
# wraps$`Fungicides and Bactericides`
# wraps$`Plant Growth Regulators`

season_cor_plots$Insecticides$use_per_land_area / 
  season_cor_plots$Herbicides$use_per_land_area / 
  season_cor_plots$`Fungicides and Bactericides`$use_per_land_area / 
  season_cor_plots$`Plant Growth Regulators`$use_per_land_area +
  plot_layout(axis_titles = 'collect')

tibble(pesticide = c(rep("insects", 3), rep("herbs", 3), rep("fungi", 3), rep("plant_regs", 3)),
       season = rep(c("spr", "sum", "aut"), 4)) %>%
  mutate(p_value = c(season_cor_tests$Insecticides %>% unlist(),
                     season_cor_tests$Herbicides %>% unlist(),
                     season_cor_tests$`Fungicides and Bactericides` %>% unlist(),
                     season_cor_tests$`Plant Growth Regulators` %>% unlist()
                     )) %>%
  mutate(p_adjust = p.adjust(p_value, method = "BH"))


# 3. As Insecticides correlate significantly, look closer
# wraps$`Insecticides – Chlorinated Hydrocarbons`
# wraps$`Insecticides – Organo-phosphates`
# wraps$`Insecticides – Carbamates`
# wraps$`Insecticides – Pyrethroids`
# wraps$`Insecticides - Biopesticides`
# wraps$`Insecticides - nes`

season_cor_plots$`Insecticides – Chlorinated Hydrocarbons`$use_per_land_area / 
  season_cor_plots$`Insecticides – Organo-phosphates`$use_per_land_area /
  season_cor_plots$`Insecticides – Carbamates`$use_per_land_area / 
  season_cor_plots$`Insecticides – Pyrethroids`$use_per_land_area / 
  season_cor_plots$`Insecticides - Biopesticides`$use_per_land_area /
  season_cor_plots$`Insecticides - nes` +
  plot_layout(axis_titles = "collect")

tibble(pesticide = c(rep("chlori", 3), rep("organo", 3), rep("carba", 3), rep("pyre", 3), rep("bio", 3), rep("nes", 3)),
       season = rep(c("spr", "sum", "aut"), 6)) %>%
  mutate(p_value = c(season_cor_tests$`Insecticides – Chlorinated Hydrocarbons` %>% unlist(),
                     season_cor_tests$`Insecticides – Organo-phosphates` %>% unlist(),
                     season_cor_tests$`Insecticides – Carbamates` %>% unlist(),
                     season_cor_tests$`Insecticides – Pyrethroids` %>% unlist(),
                     season_cor_tests$`Insecticides - Biopesticides` %>% unlist(),
                     season_cor_tests$`Insecticides - nes` %>% unlist())) %>%
  mutate(p_adjust = p.adjust(p_value, method = "BH"))

# 4. Bonus
season_cor_plots$`Herbicides – Triazines`
season_cor_plots$`Fung & Bact – Benzimidazoles`



#


pest_cor_c_s_tibble <- sulf_meta %>%
  filter(sulf_metabolism) %>%
  group_by(Country, Season) %>%
  mutate(c_s_mean_sulf_tpm = mean(sulf_tpm)) %>%
  filter(c_s_mean_sulf_tpm > 0) %>%
  ungroup() %>%
  select(Country, Season, c_s_mean_sulf_tpm) %>%
  distinct() %>%
  arrange(desc(c_s_mean_sulf_tpm)) %>%
  full_join(., pest_cor_tibble, by = "Country")

c_s_cor_tests <- list()
c_s_cor_plots <- list()
for (season in c("spr", "sum", "aut")) {
  for (year in c("use_per_cap_2019", "use_per_cap_2020")) {
    bla <- pest_cor_c_s_tibble %>%
      filter(Season == season) %>%
      select(Country, Season, c_s_mean_sulf_tpm, all_of(year))
    c_s_cor_tests[[year]][[season]] <- cor.test(bla$c_s_mean_sulf_tpm, bla[[year]])
    
    c_s_cor_plots[[year]][[season]] <- bla %>%
      ggplot(aes(x = .data[[year]], y = c_s_mean_sulf_tpm)) +
      geom_point() +
      geom_smooth(method = "glm", formula = y ~ x) +
      stat_cor(method = "spearman", label.x = 0, label.y = max(bla$c_s_mean_sulf_tpm)*1.1) +
      ggtitle(paste0(year, " - ", season))
    
  }
}
c_s_wrap <- wrap_plots(wrap_plots(c_s_cor_plots$use_per_cap_2019) /wrap_plots(c_s_cor_plots$use_per_cap_2020))



pest_use_radius <- FAOSTAT_added_data %>%
  filter(Year == 2019) %>%
  select(Country, Item, `Use per area of cropland`) %>%
  distinct() %>%
  left_join(., cropland_fraction, by = "Country") %>%
  select(-Latitude, -Longitude) %>%
  mutate(pest_use_2km = `Use per area of cropland` * cropland_fraction_2k_radius)

radius_plots <- list()
radius_p_values <- list()
for (item in unique(pest_use_radius$Item)) {
  for (season in c("spr", "sum", "aut")) {
    radius_tibble <- pest_use_radius %>%
      filter(Item == item) %>%
      # left_join(., sulf_stats$Country, by = "Country")
      left_join(., pest_cor_c_s_tibble, by = "Country") %>%
      filter(Season  == season)

    
    radius_p_values[[item]][[season]] <- cor.test(radius_tibble$pest_use_2km, radius_tibble$c_s_mean_sulf_tpm, method = "spearman")
    
    radius_plots[[item]][[season]] <- radius_tibble %>%
      ggplot(aes(x = pest_use_2km, y = c_s_mean_sulf_tpm)) +
      geom_point() +
      geom_smooth(method = "glm", formula = y ~ x) +
      stat_cor(method = "spearman", label.x = 0, label.y = max(radius_tibble$c_s_mean_sulf_tpm)*1.1) +
      ggtitle(paste0(item, " - ", season))
  }
}

pest_use_radius %>%
  filter(Item == "Pesticides (total)") %>%
  left_join(., sulf_stats$Country, by = "Country") %>%
  ggplot(aes(x = pest_use_2km, y = mean_sulf_tpm)) +
  geom_point() +
  geom_smooth(method = "glm", formula = y ~ x) +
  stat_cor(method = "spearman") +
  ggtitle("all pests all seasons")
wrap_plots(radius_plots$`Pesticides (total)`) + plot_layout(axis_titles = "collect")
wrap_plots(radius_plots$Insecticides) + plot_layout(axis_titles = "collect")
wrap_plots(radius_plots$Herbicides) + plot_layout(axis_titles = "collect")
wrap_plots(radius_plots$`Fungicides and Bactericides`) + plot_layout(axis_titles = "collect")
wrap_plots(radius_plots$`Plant Growth Regulators`) + plot_layout(axis_titles = "collect")

sulf_meta %>%
  left_join(., classification, by = "contig") %>%
  select(contig, Core, sulf_metabolism) %>%
  distinct() %>%
  group_by(sulf_metabolism, Core) %>%
  summarise(counts = n()) %>%
  group_by(Core) %>%
  mutate(prev = counts / sum(counts)) %>%
  arrange(Core)

## Save files
system("mkdir -p output/R/gene_content/sulfur")

write_delim(sulf_positive_hives, "output/R/gene_content/sulfur/sulfur_positive_hives.tsv", delim = "\t ")
for (meta in names(sulf_stats)) {
  write_delim(sulf_stats[[meta]], paste0("output/R/gene_content/sulfur/sulfur_stats.", meta, ".tsv"), delim = "\t")
  write_delim(tpm_KW[[meta]], paste0("output/R/gene_content/sulfur/tpm_krusk.", meta, ".tsv"), delim = "\t")
  write_delim(genomes_KW[[meta]], paste0("output/R/gene_content/sulfur/genome_count_krusk.", meta, ".tsv"), delim = "\t")
}

for (meta in names(sulf_tpm_plots)) {
  ggsave(paste0("output/R/gene_content/sulfur/tpm.", meta, ".pdf"),
         sulf_tpm_plots[[meta]], width = 7, height = 7)
  ggsave(paste0("output/R/gene_content/sulfur/genome_count.", meta, ".pdf"),
         sulf_genome_prop_plots[[meta]], width = 7, height = 7)
  for (set in names(nosema_sulf_cor_plot)) {
    ggsave(paste0("output/R/gene_content/sulfur/nosema_mapped_reads_prop.", meta, ".", set, ".pdf"),
           nosema_mapped_prop_plots[[meta]][[set]], width = 7, height = 7)
    
  }
}

for (set in names(nosema_sulf_cor_plot)) {
  write_delim(nosema_sulf_cor_tibble[[set]], paste0("output/R/gene_content/sulfur/nosema_sulfur_correlation.", set, ".tsv"), delim = "\t")
  ggsave(paste0("output/R/gene_content/sulfur/nosema_sulfur_correlation.", set, ".pdf"),
         nosema_sulf_cor_plot[[set]], width = 6, height = 6)
}


ggsave("output/R/gene_content/sulfur/pest_cor_plot.pdf", pest_cor_plot, width = 8, height = 6)
ggsave("output/R/gene_content/sulfur/pest_c_s_cor_plot.pdf", c_s_wrap, width = 12, height = 8)

