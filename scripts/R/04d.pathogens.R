library(readxl)
library(lme4)
library(lmerTest)
library(ggpubr)
library(tidyverse)


metadata <- readRDS("output/R/R_variables/metadata.RDS") %>% 
  tibble() %>%
  mutate(Hive_ID = as.character(Hive_ID))


gene_tpm <- read.delim("output/R/gene_content/landuse/gene_tpm.tsv") %>% tibble()

pathogen_data <- read_excel("data/GlobalBGOOD_WP1_Tier1_Scien.xlsx", skip = 1) %>%
  rename(BGOOD_sample_code = Sample_ID) %>%
  rename(`Cat DWV A` = Cat....24,
         `Cat DWV B` = Cat....26,
         `Cat ABPV` = Cat....28,
         `Cat CBPV` = Cat....30,
         `Cat BQCV` = Cat....32,
         `Cat SBV` = Cat....34,
         `Cat EFB` = Cat....36,
         `Cat AFB` = Cat....38) %>%
  mutate(across(starts_with("Cat"), ~ na_if(.x, "-")),
         across(starts_with("Cat"), ~ factor(.x, levels = c("L", "M", "H"))))

pathogens_Cts <- c("DWV A", "DWV B", "ABPV", "CBPV", "BQCV", "SBV", "EFB",
               "AFB", "AFB (cfu)", "N. apis", "N. ceranae", "N. spores")

pathogens_cats <- c("Cat DWV A", "Cat DWV B", "Cat ABPV", "Cat CBPV", "Cat BQCV",
               "Cat SBV", "Cat EFB", "Cat AFB")

meta_with_pathogens <- pathogen_data %>% 
  filter(Years == 2020) %>%
  select(BGOOD_sample_code, all_of(c(pathogens_Cts, pathogens_cats))) %>%
  left_join(metadata[c("Sample_ID", "BGOOD_sample_code", "Bee_pool", "Country", "Hive_ID", "Season", "Gut_part")], ., by = "BGOOD_sample_code") %>%
  distinct() %>%
  mutate(across(all_of(pathogens_Cts), ~ str_replace_all(., "negative", "40"))) %>%
  mutate(across(all_of(pathogens_Cts), ~ suppressWarnings(as.numeric(.))))

mean_gene_tpm_and_pathogens <- gene_tpm %>% 
  left_join(., meta_with_pathogens, by = "Sample_ID") %>%
  # group_by(gene, Bee_pool) %>%
  # mutate(mean_tpm = mean(tpm), .before = BGOOD_sample_code) %>% 
  # ungroup() %>%
  select(gene, Sample_ID, Country, Hive_ID, Season, Gut_part, tpm, all_of(c(pathogens_Cts, pathogens_cats))) %>%
  distinct()

informative_numbers <- tibble()
for (pathogen in c(pathogens_Cts, pathogens_cats)) {
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

gene_pathogen_cor_tibble <- list()
model <- list()
coeffs <- tibble()

pathogens_of_interest <- c("BQCV", "SBV", "DWV B")
# for (goi in unique(gene_tpm$gene)) {
for (goi in "phosphoadenosine phosphosulfate reductase") {
  for (poi in pathogens_of_interest) {
    # goi = "phosphoadenosine phosphosulfate reductase"
    # poi = "BQCV"
    
    gene_pathogen_cor_tibble[[goi]][[poi]] <- mean_gene_tpm_and_pathogens %>%
      filter(gene == goi) %>%
      rename(pathogen_ct = all_of(poi)) %>% 
      filter(!is.na(pathogen_ct)) %>%
      select(Sample_ID, Country, Hive_ID, Season, Gut_part, tpm, pathogen_ct) %>%
      # mutate(logged_ct = log2(pathogen_ct)) %>%
      mutate(log_tpm = log10(tpm)) %>%
      filter(!is.infinite(log_tpm))
    
    model[[goi]][[poi]] <- lmer(log_tpm ~ pathogen_ct + Season + Gut_part + 
           (1 | Hive_ID ),
         data = gene_pathogen_cor_tibble[[goi]][[poi]])
    
    coeffs <- summary(model[[goi]][[poi]])$coefficients %>%
      as.data.frame() %>%
      rownames_to_column("metric") %>%
      tibble() %>%
      mutate(gene = goi, 
             pathogen = poi,
             .before = metric) %>%
      mutate(singular = ifelse(isSingular(model[[goi]][[poi]]), TRUE, FALSE)) %>%
      rbind(coeffs)
  }
}

summary(model$`phosphoadenosine phosphosulfate reductase`$BQCV)
summary(model$`phosphoadenosine phosphosulfate reductase`$SBV)
summary(model$`phosphoadenosine phosphosulfate reductase`$`DWV B`)

plot(model$`phosphoadenosine phosphosulfate reductase`$BQCV, which = 1)
qqnorm(resid(model$`phosphoadenosine phosphosulfate reductase`$BQCV))


slopes <- coeffs %>%
  filter(metric == "pathogen_ct") %>%
  mutate(p_adjusted = p.adjust(`Pr(>|t|)`, method = "BH")) %>%
  mutate(p_adj_sig = case_when(p_adjusted <= 0.001 ~ "***",
                               p_adjusted <= 0.01 ~ "**",
                               p_adjusted <= 0.05 ~ "*",
                               p_adjusted <= 0.075 ~ ".",
                               .default = "n.s."
                               ))



goi = "phosphoadenosine phosphosulfate reductase"
poi = "Cat BQCV"

gene_pathogen_cor_tibble[[goi]][[poi]] <- mean_gene_tpm_and_pathogens %>%
  filter(gene == goi) %>%
  rename(pathogen_ct = all_of(poi)) %>% 
  filter(!is.na(pathogen_ct)) %>%
  select(Bee_pool, Country, Hive_ID, Season, mean_tpm, pathogen_ct) %>%
  # mutate(logged_ct = log2(pathogen_ct)) %>%
  mutate(log_mean_tpm = log10(mean_tpm)) %>%
  filter(!is.infinite(log_mean_tpm))


gene_pathogen_cor_tibble[[goi]][[poi]] %>%
  ggplot(aes(x = pathogen_ct, y = log_mean_tpm))




#

pathogens_of_interest <- c("BQCV", "SBV", "DWV B")
# pathogens_of_interest <- c("BQCV", "SBV", "DWV B", "ABPV", "N. ceranae", "DWV A", "CBPV", "N. spores")
gene_pathogen_cor_tibble <- list()
gene_pathogen_cor_plot <- list()
gene_pathogen_cor_test <- list()
gene_pathogen_cor_p_values <- tibble()
# for (goi in unique(gene_tpm$gene)) {
for (goi in "phosphoadenosine phosphosulfate reductase") {
    for (poi in pathogens_of_interest) {
      # goi = "phosphoadenosine phosphosulfate reductase"
      # poi = "BQCV"
    gene_pathogen_cor_tibble[[goi]][[poi]] <- mean_gene_tpm_and_pathogens %>%
      filter(gene == goi) %>%
      filter(!is.na(.data[[poi]])) %>%
      select(Sample_ID, Country, Hive_ID, Season, tpm, all_of(poi)) %>%
      rename(pathogen_ct = all_of(poi)) %>%
      # mutate(logged_ct = log2(.data[[poi]])) %>%
      mutate(log_tpm = log10(tpm)) %>%
      filter(!is.infinite(log_tpm))
    
    if (nrow(gene_pathogen_cor_tibble[[goi]][[poi]]) > 2) {
      gene_pathogen_cor_plot[[goi]][[poi]] <- gene_pathogen_cor_tibble[[goi]][[poi]] %>%
        # ggplot(aes(x = .data[[poi]], y = log_mean_tpm, color = Country)) +
        # ggplot(aes(x = .data[[poi]], y = log_mean_tpm)) +
        ggplot(aes(x = pathogen_ct, y = log_tpm)) +
        # ggplot(aes(x = .data[[poi]], y = mean_tpm)) +
        geom_point() +
        geom_smooth(method = "glm", formula = y ~ x) +
        # stat_cor(method = "spearman") +
        stat_cor(method = "pearson") +
        ggtitle(paste0(goi, " - ", poi))
      gene_pathogen_cor_test[[goi]][[poi]] <- cor.test(gene_pathogen_cor_tibble[[goi]][[poi]]$pathogen_ct,  
               gene_pathogen_cor_tibble[[goi]][[poi]]$log_tpm,
               # method = "spearman")
               method = "pearson")
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
