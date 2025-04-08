library(ggpubr)
library(tidyverse)

metadata <- readRDS("output/R/R_variables/metadata.RDS") %>% 
  tibble() %>%
  mutate(Hive_ID = as.character(Hive_ID))


gene_tpm <- read.delim("output/R/gene_content/landuse/gene_tpm.tsv") %>% tibble()

pathogen_data <- read_excel("data/GlobalBGOOD_WP1_Tier1_Scien.xlsx", skip = 1) %>%
  rename(BGOOD_sample_code = Sample_ID)

pathogens <- c("DWV A", "DWV B", "ABPV", "CBPV", "BQCV", "SBV", "EFB", 
               "AFB", "AFB (cfu)", "N. apis", "N. ceranae", "N. spores")

meta_with_pathogens <- pathogen_data %>% 
  filter(Years == 2020) %>%
  select(BGOOD_sample_code, all_of(pathogens)) %>%
  left_join(metadata[c("Sample_ID", "BGOOD_sample_code", "Bee_pool", "Country", "Hive_ID", "Season", "Gut_part")], ., by = "BGOOD_sample_code") %>%
  distinct() %>%
  mutate(across(all_of(pathogens), ~ str_replace_all(., "negative", "0"))) %>%
  mutate(across(all_of(pathogens), ~ suppressWarnings(as.numeric(.))))

mean_gene_tpm_and_pathogens <- gene_tpm %>% 
  left_join(., meta_with_pathogens, by = "Sample_ID") %>%
  group_by(gene, Bee_pool) %>%
  mutate(mean_tpm = mean(tpm), .before = BGOOD_sample_code) %>% 
  ungroup() %>%
  select(gene, Bee_pool, Country, Hive_ID, Season, mean_tpm, all_of(pathogens)) %>%
  distinct()

informative_numbers <- tibble()
for (pathogen in pathogens) {
  NAs <- mean_gene_tpm_and_pathogens %>%
    filter(is.na(.data[[pathogen]])) %>%
    nrow()
  zeros <- mean_gene_tpm_and_pathogens %>%
    filter(.data[[pathogen]] == 0) %>%
    nrow()
  informative_numbers <- tibble(pat = pathogen, 
         NA_prop = NAs / nrow(mean_gene_tpm_and_pathogens), 
         zero_prop = zeros / nrow(mean_gene_tpm_and_pathogens),
         informative_numbers_prop = (nrow(mean_gene_tpm_and_pathogens) - NAs - zeros) / nrow(mean_gene_tpm_and_pathogens)) %>%
    rbind(informative_numbers, .)
}

informative_numbers %>%
  arrange(desc(informative_numbers_prop))

pathogens_of_interest <- c("BQCV", "SBV", "DWV B")
gene_pathogen_cor_tibble <- list()
gene_pathogen_cor_plot <- list()
gene_pathogen_cor_test <- list()
gene_pathogen_cor_p_values <- tibble()
for (goi in unique(gene_tpm$gene)) {
  for (poi in pathogens_of_interest) {
    gene_pathogen_cor_tibble[[goi]][[poi]] <- mean_gene_tpm_and_pathogens %>%
      filter(gene == goi) %>%
      filter(!is.na(.data[[poi]])) %>%
      select(Bee_pool, Country, Hive_ID, Season, mean_tpm, all_of(poi)) %>%
      mutate(log_mean_tpm = log10(mean_tpm)) %>%
      filter(!is.infinite(log_mean_tpm))
    
    if (nrow(gene_pathogen_cor_tibble[[goi]][[poi]]) > 2) {
      gene_pathogen_cor_plot[[goi]][[poi]] <- gene_pathogen_cor_tibble[[goi]][[poi]] %>%
        # ggplot(aes(x = .data[[poi]], y = log_mean_tpm, color = Country)) +
        ggplot(aes(x = .data[[poi]], y = log_mean_tpm)) +
        geom_point() +
        geom_smooth(method = "glm", formula = y ~ x) +
        stat_cor(method = "spearman") +
        ggtitle(paste0(goi, " - ", poi))
      gene_pathogen_cor_test[[goi]][[poi]] <- cor.test(gene_pathogen_cor_tibble[[goi]][[poi]][[poi]],  
               gene_pathogen_cor_tibble[[goi]][[poi]]$mean_tpm,
               method = "spearman")
      gene_pathogen_cor_p_values <- tibble(
        test = paste0(goi, " - ", poi),
        R = gene_pathogen_cor_test[[goi]][[poi]]$estimate,
        p_value = gene_pathogen_cor_test[[goi]][[poi]]$p.value) %>%
        rbind(gene_pathogen_cor_p_values, .)
    }
  }
}

gene_pathogen_cor_p_values %>%
  mutate(adjusted_p = p.adjust(p_value, method = "BH")) %>%
  arrange(adjusted_p)


gene_pathogen_cor_p_values %>%
  filter(str_detect(test, "phosphoadenosine phosphosulfate reductase")) %>%
  mutate(adjusted_p = p.adjust(p_value, method = "BH")) %>%
  arrange(adjusted_p)
