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

phage_tpm <- read.csv("output/R/relative_abundance/phage_tpm.csv") %>%
  tibble()

phold_predictions_with_extensions <- read.csv("output/R/gene_content/phold_predictions_with_extensions.csv") %>%
  tibble() %>%
  filter(str_starts(contig_id, "NODE"))

kegg_mapping <- read.delim("data/kegg_mapping.tsv", colClasses = "character") %>%
  tibble()

kegg_and_phold <- kegg_mapping %>%
  left_join(., phold_predictions_with_extensions[c("cds_id", "phrog", "function.", "product")], by = "cds_id")

CDSs_with_metabolism_kegg <- kegg_and_phold %>%
  filter(Pathway_category == "Metabolism" | 
           product %in% c("chitinase", "glutamine amidotransferase", 
                          "PnuC-like nicotinamide mononucleotide transport")) %>%
  filter(!product %in% c("decoy of host sigma70", "MazF-like growth inhibitor",
                         "toxin", "VFDB virulence factor protein")) %>%
  filter(!str_detect(product, "Que")) %>% # This will remove 3 genes. All of them are only present in one sample (the same one for all 3)
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
  pivot_wider(names_from = product, values_from = present, values_fill = 0)

gene_tpm <- grene_presence_on_contigs %>%
  pivot_longer(-contig, names_to = "gene", values_to = "present") %>%
  left_join(., phage_tpm, by = "contig") %>% 
  mutate(across(-c(contig, gene, present), ~ .x * present)) %>%
  select(-present) %>%
  pivot_longer(-c(contig, gene), names_to = "Sample_ID", values_to = "tpm") %>%
  group_by(gene, Sample_ID) %>%
  mutate(tpm = sum(tpm)) %>%
  ungroup() %>%
  select(-contig) %>%
  distinct()

genes_of_interest <- unique(gene_tpm$gene)
coeffs_tpm_simple <- list()
countries_left <- list()

#####
# CROPLAND
test_tibble_log_tpm <- list()
model_tpm_simple_cropland <- list()
for (goi in genes_of_interest) {
  
  test_tibble_log_tpm[[goi]] <- gene_tpm %>%
    filter(gene == goi) %>%
    left_join(., metadata[c("Sample_ID", "Country", "Hive_ID", "Season", "Gut_part")], by = "Sample_ID") %>%
    distinct() %>%
    left_join(., cropland_and_FAO, by = "Country") %>%
    mutate(Gut_part = factor(Gut_part, levels = c("rec", "ile", "mid"))) %>%
    # mutate(log_tpm = log10(tpm)) %>%
    # filter(!is.infinite(log_tpm)) %>%
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
  
  model_tpm_simple_cropland[[goi]] <- brm(tpm ~ ha_cropland_in_2k_radius + Gut_part + Season +
                                             (1 | Country / Hive_ID ), 
                                          data = test_tibble_log_tpm[[goi]],
                                          family = brmsfamily("hurdle_lognormal"),
                                          chains = 4, iter = 10000, cores = 4,
                                          seed = 1,
                                          prior = prior(normal(-3, 0.5), class = "Intercept", ub = 0),
                                          control = list(adapt_delta = 0.99)
                                          )
  
  summary(model_tpm_simple_cropland[[goi]])
  
  
}




#

x <- seq(0.001, 1, length.out = 1000)

# Calculate density using the lognormal distribution (meanlog = -1, sdlog = 0.5)

mlog <- -3.5
sdlog = 0.5
density_values <- dlnorm(x, meanlog = mlog, sdlog = sdlog)

# Put the values in a data frame
df <- data.frame(tpm = x, density = density_values)

# Plot using ggplot2
ggplot(df, aes(x = tpm, y = density)) +
  geom_line(color = "blue", size = 1) +
  labs(title = paste0("Implied Prior on TPM (Lognormal with μ = ", mlog, ", σ = ", sdlog, ")"),
       x = "TPM", y = "Density") +
  theme_minimal()

test_tibble_log_tpm[[goi]] %>% filter(tpm > 0 ) %>% ggplot(aes(x=tpm)) + geom_histogram(bins = 100)






#