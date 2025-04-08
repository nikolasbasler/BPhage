library(tidyverse)
library(ggpubr)
library(patchwork)

metadata <- readRDS("output/R/R_variables/metadata.RDS")
classification <- readRDS("output/R/R_variables/classification.RDS")
present_in_all_countries <- read_lines("data/core_contigs.txt")

# all_hosts <- read.csv("output/R/host_pies/all_hosts.all.csv") %>%
#   rename(contig = Virus, Host_genus = Genus)

phold_predictions_with_extensions <- read.csv("output/R/gene_content/phold_predictions_with_extensions.csv") %>%
  tibble() %>%
  filter(str_starts(contig_id, "NODE"))

cropland_fraction <- read.csv("data/land_cover_results.csv") %>% 
  tibble() %>%
  mutate(cropland_fraction = cropland_fraction / 100) %>%
  rename(cropland_fraction_2k_radius = cropland_fraction) %>%
  arrange(cropland_fraction_2k_radius)

phage_tpm <- read.csv("output/R/relative_abundance/phage_tpm.csv") %>%
  tibble()

genes_of_interest <- c("phosphoadenosine phosphosulfate reductase", "levanase",
                       "ribosomal protein S6 glutaminyl transferase", "glucosyltransferase",
                       "glutamine amidotransferase", "porphyrin biosynthesis", 
                       "aerobic cobaltochelatase CobT subunit", "chitinase",
                       "dTDP-4-dehydrorhamnose 3", "nicotinamide-nucleotide adenylyltransferase",
                       "PnuC-like nicotinamide mononucleotide transport",
                       "NrdD-like anaerobic ribonucleotide reductase large subunit",
                       "QueD-like  6-pyruvoyl-tetrahydropterin synthase",
                       "QueC-like queuosine biosynthesis", "QueE-like  radical SAM domain")
genes_of_particular_interest <- c("phosphoadenosine phosphosulfate reductase", "levanase")

gene_tpm <- phold_predictions_with_extensions %>%
  rename("contig" = "contig_id") %>%
  filter(function. == "moron, auxiliary metabolic gene and host takeover") %>% 
  group_by(product) %>% 
  mutate(gene_count = n()) %>%
  ungroup() %>%
  filter(gene_count >= 5) %>% # Maybe adapt this
  select(contig, product) %>% 
  left_join(., phage_tpm, by = "contig") %>%
  pivot_longer(-c(contig, product), names_to = "Sample_ID", values_to = "TPM") %>%
  distinct() %>%
  group_by(product, Sample_ID) %>%
  summarise(gene_tpm = sum(TPM), .groups = "drop") %>%
  # left_join(., metadata[c("Sample_ID", "Country", "Hive_ID", "Gut_part", "Season")], by = "Sample_ID") %>%
  left_join(., metadata[c("Sample_ID", "Country", "Gut_part")], by = "Sample_ID") %>%
  left_join(., cropland_fraction[c("Country", "cropland_fraction_2k_radius")], by = "Country")

########
#######

absolute_counts <- read.csv("output/R/absolute_counts.csv") %>%
  tibble()
kegg_mapping <- read.delim("data/kegg_mapping.tsv", colClasses = "character") %>%
  tibble()
kegg_and_phold <- kegg_mapping %>%
  left_join(., phold_predictions_with_extensions[c("cds_id", "phrog", "function.", "product")], by = "cds_id")


CDSs_with_metabolism_kegg <- kegg_and_phold %>%
  # filter(Pathway_category == "Metabolism") %>%
  # filter(Pathway_category %in% c("Metabolism", "Environmental Information Processing", NA)) %>%
  filter(Pathway_category == "Metabolism" | 
           product %in% c("chitinase", "glutamine amidotransferase", 
                          "PnuC-like nicotinamide mononucleotide transport")) %>%
  filter(!product %in% c("decoy of host sigma70", "MazF-like growth inhibitor", 
                         "toxin", "VFDB virulence factor protein")) %>%
  distinct(cds_id) %>% 
  unlist(use.names = FALSE)
genes_with_kegg <- phold_predictions_with_extensions %>% 
  filter(cds_id %in% CDSs_with_metabolism_kegg) %>% 
  distinct(product) %>%
  unlist(use.names = FALSE)

grene_presence_on_contigs <- phold_predictions_with_extensions %>% 
  filter(product %in% genes_with_kegg) %>%
  rename(contig = contig_id) %>%
  select(contig, product) %>% 
  mutate(present = 1) %>%
  distinct() %>% 
  # mutate(product = str_replace_all(product, "-", "_"),
  #        product = str_replace_all(product, " ", "_")) %>%
  pivot_wider(names_from = product, values_from = present, values_fill = 0)


gene_load_long <- grene_presence_on_contigs %>%
  filter(contig %in% absolute_counts$contig) %>%
  pivot_longer(-contig, names_to = "gene", values_to = "present") %>%
  left_join(., absolute_counts, by = "contig") %>%
  mutate(across(-c(contig, gene, present), ~ .x * present)) %>%
  select(-present) %>%
  pivot_longer(-c(contig, gene), names_to = "Sample_ID", values_to = "viral_load") %>%
  group_by(gene, Sample_ID) %>%
  mutate(viral_load = sum(viral_load)) %>%
  ungroup() %>%
  select(-contig) %>%
  left_join(., metadata[c("Sample_ID", "Country", "Gut_part")], by = "Sample_ID") %>%
  mutate(log_load = log(viral_load)) %>%
  filter(!is.infinite(log_load)) %>%
  left_join(., cropland_fraction[c("Country", "cropland_fraction_2k_radius")], by = "Country")
  

#######
#######

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
         Year == 2019) %>%
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

gut_part_gene_cor_plots <- list()
raw_p_values <- list()
# for (g_part in unique(gene_tpm$Gut_part)) {
for (g_part in "rec") {
  # for (prod in unique(gene_tpm$product)) {
  for (prod in genes_of_interest) {
    for (item in c("Pesticides (total)", "Insecticides", "Herbicides", "Fungicides and Bactericides", "Plant Growth Regulators")) {
      for (element in c("Use per area of cropland", "Use per capita", "Use per value of agricultural production", "use_per_active_cropland", "use_per_land_area", "est_use_in_2k_radius")) {
        # plot_tible <- gene_tpm %>%
        #   filter(product == prod,
        #          Gut_part == g_part) %>% 
        #   # mutate(log_gene_tpm = log(gene_tpm)) %>%
        #   group_by(Country) %>%
        #   mutate(mean_gene_tpm = mean(gene_tpm)) %>%
        #   # filter(!is.infinite(log_gene_tpm)) %>%
        #   # mutate(mean_gene_tpm = mean(log_gene_tpm)) %>%
        #   ungroup() %>%
        #   select(Country, product, mean_gene_tpm, cropland_fraction_2k_radius) %>%
        #   distinct() %>%
        #   left_join(., FAOSTAT_added_data, by = "Country") %>%
        #   filter(Item == item) %>%
        #   select(Country, product, mean_gene_tpm, all_of(element), cropland_fraction_2k_radius)
        plot_tible <- gene_load_long %>%
          filter(gene == prod) %>%
          group_by(Country) %>%
          mutate(mean_log_load = mean(log_load)) %>%
          ungroup() %>%
          select(Country, gene, mean_log_load, cropland_fraction_2k_radius) %>%
          distinct() %>%
          left_join(., FAOSTAT_added_data, by = "Country") %>%
          filter(Item == item)
        
        if (nrow(plot_tible) < 2) {
          next
        }
        gut_part_gene_cor_plots[[g_part]][[prod]] <- plot_tible %>%
          # ggplot(aes(x = .data[[element]], y = mean_gene_tpm)) +
          ggplot(aes(x = .data[[element]], y = mean_log_load)) +
          geom_point() +
          geom_smooth(method = "glm", formula = y ~ x) +
          stat_cor(method = "spearman") +
          ggtitle(paste0(g_part, " - ", prod, " - ", item, " - ", element))

        # raw_p_values[[paste(g_part, prod, item, element, sep = "_")]] <- cor.test(plot_tible[[element]], plot_tible$mean_gene_tpm, method = "spearman")$p.value 
        raw_p_values[[paste(g_part, prod, item, element, sep = "_")]] <- cor.test(plot_tible[[element]], plot_tible$mean_log_load, method = "spearman")$p.value 
      }
    }
    # raw_p_values[[paste(g_part, prod, "cropland_fraction_2k_radius", "NA", sep = "_")]] <- cor.test(plot_tible$cropland_fraction_2k_radius, plot_tible$mean_gene_tpm, method = "spearman")$p.value
  }
}

gut_part_gene_cor_plots$mid$`phosphoadenosine phosphosulfate reductase`
gut_part_gene_cor_plots$ile$`phosphoadenosine phosphosulfate reductase`
gut_part_gene_cor_plots$rec$`phosphoadenosine phosphosulfate reductase`

as_tibble(raw_p_values) %>% 
  rownames_to_column("row") %>% 
  pivot_longer(-row, names_to = "test", values_to = "p_value") %>%
  select(-row) %>%
  separate_wider_delim(test, "_", names = c("Gut_part", "product", "pesticide","metric"), too_many = "merge") %>% 
  filter(p_value <= 0.05,
         metric %in% c("Use per area of cropland", "use_per_land_area", "est_use_in_2k_radius", "fraction_2k_radius_NA")) %>% View()
