library(lme4)
library(lmerTest)
# library(glmmTMB)
library(ggpubr)
library(DHARMa)
library(tidyverse)

cropland_fraction <- read.csv("data/land_cover_results.csv") %>% 
  tibble() %>%
  mutate(cropland_fraction = cropland_fraction / 100) %>%
  rename(cropland_fraction_2k_radius = cropland_fraction) %>%
  arrange(cropland_fraction_2k_radius)
# make FAOSTAT_added_data:
source("scripts/R/helpers/FAOstat_table.R")

metadata <- readRDS("output/R/R_variables/metadata.RDS") %>%
  mutate(Hive_ID = as.character(Hive_ID))
classification <- readRDS("output/R/R_variables/classification.RDS")

phage_tpm <- read.csv("output/R/relative_abundance/phage_tpm.csv") %>%
  tibble()

absolute_counts <- read.csv("output/R/absolute_counts.csv") %>%
  tibble()

phold_predictions_with_extensions <- read.csv("output/R/gene_content/phold_predictions_with_extensions.csv") %>%
  tibble() %>%
  filter(str_starts(contig_id, "NODE"))

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
  select(-contig)

gene_tpm_long <- list()
gene_tpm_long$Sample_ID <- grene_presence_on_contigs %>%
  pivot_longer(-contig, names_to = "gene", values_to = "present") %>%
  left_join(., phage_tpm, by = "contig") %>% 
  mutate(across(-c(contig, gene, present), ~ .x * present)) %>%
  select(-present) %>%
  pivot_longer(-c(contig, gene), names_to = "Sample_ID", values_to = "tpm") %>%
  group_by(gene, Sample_ID) %>%
  mutate(tpm = sum(tpm)) %>%
  ungroup() %>%
  select(-contig)
# 
# gene_tpm_long$Bee_pool <- gene_tpm_long$Sample_ID %>%
#   left_join(., metadata[c("Sample_ID", "Bee_pool")], by = "Sample_ID") %>%
#   group_by(gene, Bee_pool) %>%
#   mutate(mean_tpm = mean(tpm)) %>%
#   ungroup() %>%
#   select(gene, Bee_pool, mean_tpm) %>%
#   distinct()
# 
# gene_tpm_long$Hive <- gene_tpm_long$Sample_ID %>%
#   left_join(., metadata[c("Sample_ID", "Hive_ID")], by = "Sample_ID") %>%
#   group_by(gene, Hive_ID) %>%
#   mutate(mean_tpm = mean(tpm)) %>%
#   ungroup() %>%
#   select(gene, Hive_ID, mean_tpm) %>%
#   distinct()
# 
# gene_tpm_long$Hive_Gut <- gene_tpm_long$Sample_ID %>%
#   left_join(., metadata[c("Sample_ID", "Country", "Hive_ID", "Gut_part")], by = "Sample_ID") %>%
#   group_by(gene, Hive_ID, Gut_part) %>%
#   mutate(mean_tpm = mean(tpm)) %>%
#   ungroup() %>%
#   mutate(Hive_Gut = paste(Hive_ID, Gut_part, sep = "_")) %>%
#   select(gene, Hive_Gut, mean_tpm, Country, Hive_ID, Gut_part) %>%
#   distinct()


genes_of_interest <- unique(gene_tpm_long$Sample_ID$gene)

# test_tibble_tpm <- list()
# model_tpm <- list()

test_tibble_tpm_log <- list()
model_tpm_log <- list()

test_tibble_abs <- list()
model_abs <- list()

# model_logit_with_zi <- list()
# model_logit <- list()
for (goi in genes_of_interest) {

  # test_tibble <- gene_tpm_long$Bee_pool %>%
  #     filter(gene == goi) %>%
  #   # filter(mean_tpm != 0)
  #   left_join(., metadata[c("Bee_pool", "Country", "Hive_ID", "Season")], by = "Bee_pool") %>%
  #   distinct() %>%
  #   # mutate(Hive_ID = as.character(Hive_ID)) %>% # CHECK IF THIS IS REALLY TREATED AS CATEGORICAL!
  #   mutate(Hive_ID = paste0("Hive_", Hive_ID)) %>%
  #   left_join(., cropland_fraction[c("Country", "cropland_fraction_2k_radius")], by = "Country")
  # 
  # model <- lmer(mean_tpm ~ cropland_fraction_2k_radius + Country + Season + (1 | Hive_ID), data = test_tibble)
  # summary(model)
  
  ### TPM
  # test_tibble_tpm[[goi]] <- gene_tpm_long$Sample_ID %>%
  #   filter(gene == goi) %>%
  #   left_join(., metadata[c("Sample_ID", "Country", "Hive_ID", "Season", "Gut_part")], by = "Sample_ID") %>%
  #   distinct() %>%
  #   # mutate(Hive_ID = paste0("Hive_", Hive_ID)) %>%
  #   left_join(., cropland_fraction[c("Country", "cropland_fraction_2k_radius")], by = "Country") %>%
  #   mutate(Gut_part = factor(Gut_part, levels = c("rec", "ile", "mid"))) # %>%
  #   # mutate(stand_cent_crop = (cropland_fraction_2k_radius - mean(cropland_fraction_2k_radius) / sd(cropland_fraction_2k_radius)))
  # 
  # model_tpm[[goi]] <- lmer(tpm ~ cropland_fraction_2k_radius * Gut_part + Season +
  #                  (1 | Hive_ID), data = test_tibble_tpm[[goi]])
  # 
  # 
  test_tibble_tpm_log[[goi]] <- gene_tpm_long$Sample_ID %>%
    filter(gene == goi) %>%
    left_join(., metadata[c("Sample_ID", "Country", "Hive_ID", "Season", "Gut_part")], by = "Sample_ID") %>% 
    distinct() %>%
    left_join(., cropland_fraction[c("Country", "cropland_fraction_2k_radius")], by = "Country") %>%
    mutate(Gut_part = factor(Gut_part, levels = c("rec", "ile", "mid"))) %>%
    mutate(log_tpm = log(tpm)) %>%
    filter(!is.infinite(log_tpm))
  
  hive_ids <- test_tibble_tpm_log[[goi]] %>% 
    distinct(Hive_ID) %>%
    nrow()
  countries_left <- test_tibble_tpm_log[[goi]] %>%
    distinct(Country) %>%
    nrow()
  # if (hive_ids < nrow(test_tibble_tpm_log[[goi]])) {
  if (countries_left > 7 ) {
    model_tpm_log[[goi]] <- lmer(log_tpm ~ cropland_fraction_2k_radius + Gut_part + Season +
                                   (1 | Hive_ID), data = test_tibble_tpm_log[[goi]])
    # model_tpm_log[[goi]] <- lmer(log_tpm ~ cropland_fraction_2k_radius * Gut_part + 
    #        cropland_fraction_2k_radius * Season + (1 | Hive_ID), 
    #      data = test_tibble_tpm_log[[goi]])
    
  } else {
    print(paste0(goi, " - model for log tpm not possible. Hive_IDs left: ", hive_ids, ". Samples left: ", nrow(test_tibble_tpm_log[[goi]]), "."))
  }
  
  # model_tpm_log[[goi]] <- lmer(log_tpm ~ cropland_fraction_2k_radius * Gut_part + Season +
  #                                (1 | Hive_ID), data = test_tibble_tpm_log[[goi]])

  ## Absolute loads

  test_tibble_abs[[goi]] <- gene_load_long %>%
    filter(gene == goi) %>%
    left_join(., metadata[c("Sample_ID", "Country", "Hive_ID", "Season")], by = "Sample_ID") %>%
    distinct() %>%
    # mutate(Hive_ID = paste0("Hive_", Hive_ID)) %>%
    left_join(., cropland_fraction[c("Country", "cropland_fraction_2k_radius")], by = "Country") %>%
    mutate(log_load = log(viral_load)) %>%
    filter(!is.infinite(log_load))

  hive_ids <- test_tibble_abs[[goi]] %>%
    distinct(Hive_ID) %>%
    nrow()
  countries_left <- test_tibble_abs[[goi]] %>%
    distinct(Country) %>%
    nrow()
  # if (hive_ids < nrow(test_tibble_tpm_log[[goi]])) {
  if (countries_left > 7 ) {
    model_abs[[goi]] <- lmer(log_load ~ cropland_fraction_2k_radius + Season +
                               (1 | Hive_ID), data = test_tibble_abs[[goi]])
    # model_abs[[goi]] <- lmer(log_load ~ cropland_fraction_2k_radius +
    #                                cropland_fraction_2k_radius * Season + (1 | Hive_ID),
    #                              data = test_tibble_abs[[goi]])

  } else {
    print(paste0(goi, " - model for absolute counts not possible. Hive_IDs left: ", hive_ids, ". Samples left: ", nrow(test_tibble_abs[[goi]]), "."))
  }

  
  # glmmTMB(viral_load ~ cropland_fraction_2k_radius + Season + (1 | Hive_ID),
  #         data = test_tibble,
  #         family = gaussian(link = "log"),
  #         ziformula = ~1)
  
  
  
  ### TPM logit
  # test_tibble_logit <- test_tibble %>%
  #   mutate(logit_tpm = log(tpm/(1-tpm))) %>%
  #   filter(!is.infinite(logit_tpm))
  # if (nrow(test_tibble_logit) < 2) {
  #   next
  # }
  # # model_logit[[goi]] <- lmer(logit_tpm ~ cropland_fraction_2k_radius * Gut_part + Season +
  # #                              (1 | Hive_ID), data = test_tibble_logit)
  # # 
  # ### TPM logit zero inflation
  # test_tibble_for_glmm <- test_tibble %>%
  #   filter(tpm != 1)
  # model_logit_with_zi[[goi]] <- glmmTMB(tpm ~ cropland_fraction_2k_radius * Gut_part + Season + (1 | Hive_ID),
  #                         data = test_tibble_for_glmm,
  #                         family = beta_family(link = "logit"),
  #                         ziformula = ~1)
  # model_logit_with_zi[[goi]] <- glmmTMB(tpm ~ cropland_fraction_2k_radius + Season + (1 | Hive_ID),
  #                                       data = test_tibble_for_glmm,
  #                                       family = beta_family(link = "logit"),
  #                                       ziformula = ~1)
  
  



  
  # summary(model[[goi]])
  
  # model2[[goi]] <- lmer(tpm ~ stand_cent_crop * Gut_part + Season +
  #                        (1 | Hive_ID), data = test_tibble)
  # summary(model[[goi]])
  
  # Model testing
  # x=summary(model[[goi]])
  # x$residuals
  
  
  # lmer(tpm ~ cropland_fraction_2k_radius +
  #        (1 | Hive_ID) + (1 | Season) + (1 | Gut_part), data = test_tibble)
 



  

  
  ### 
  # 
  # test_tibble <- gene_tpm_long$Hive %>%
  #   filter(gene == goi) %>%
  #   left_join(., metadata[c("Hive_ID", "Country")], by = "Hive_ID") %>% 
  #   distinct() %>%
  #   mutate(Hive_ID = paste0("Hive_", Hive_ID)) %>%
  #   left_join(., cropland_fraction[c("Country", "cropland_fraction_2k_radius")], by = "Country")
  # 
  # model <- lmer(mean_tpm ~ cropland_fraction_2k_radius + (1 | Country), data = test_tibble)
  # summary(model)
  # 
  # ###
  # 
  # test_tibble <- gene_tpm_long$Hive_Gut %>%
  #   filter(gene == goi) %>%
  #   left_join(., cropland_fraction, by = "Country") %>%
  #   select(Hive_Gut, mean_tpm, Hive_ID, Gut_part, cropland_fraction_2k_radius)
  # 
  # model <- lmer(mean_tpm ~ cropland_fraction_2k_radius + Gut_part + (1 | Hive_ID), data = test_tibble)
  # summary(model)
  # 
  # model <- lmer(mean_tpm ~ cropland_fraction_2k_radius * Gut_part + (1 | Hive_ID), data = test_tibble)
  # summary(model)

  ###

  # 
  # 
  
  ###
  # 
  # lmer(mean_tpm ~ cropland_fraction_2k_radius + Country + Season + (1 | Hive_ID), data = test_tibble)
  # 
  # lmer(mean_tpm ~ cropland_fraction_2k_radius + Country + Season + 
  #        (1 + cropland_fraction_2k_radius | Hive_ID), 
  #      data = test_tibble, REML = FALSE)
  # 
  # lmer(mean_tpm ~ cropland_fraction_2k_radius + Season + (1 | Hive_ID) + (1 | Country), data = test_tibble)
  # 
  # model <- lmer(mean_tpm ~ cropland_fraction_2k_radius + Season + (1 | Hive_ID), data = test_tibble)
  # summary(model)
  # isSingular(model)
  # 
  # model_simple <- lm(mean_tpm ~ cropland_fraction_2k_radius + Season, data = test_tibble)
  # summary(model_simple)
  # 
  # model2 <- lmer(tpm ~ cropland_fraction_2k_radius + Season + Gut_part + (1 | Hive_ID/Season), data = test_tibble)
  # summary(model2)
}

coeffs_cropland_tpm_log <- tibble()
for (mod in names(model_tpm_log)) {
  print(paste0("Gene: ", mod, " tpm"))
  sing <- FALSE
  if (isSingular(model_tpm_log[[mod]])) {
    print("MODEL IS SINGULAR!")
    sing <- TRUE
  }
  print(summary(model_tpm_log[[mod]]))
  coeffs_cropland_tpm_log <- summary(model_tpm_log[[mod]])$coefficients %>%
    as.data.frame() %>%
    rownames_to_column("metric") %>%
    tibble() %>%
    mutate(gene = mod, .before = metric) %>%
    mutate(singular = sing) %>%
    rbind(coeffs_cropland_tpm_log)
  
  print("")
  print("----------------------------------------------")
  print("----------------------------------------------")
  print("")
}

coeffs_cropland_abs <- tibble()
for (mod in names(model_abs)) {
  print(paste0("Gene: ", mod, " tpm"))
  sing <- FALSE
  if (isSingular(model_abs[[mod]])) {
    print("MODEL IS SINGULAR!")
    sing <- TRUE
  }
  print(summary(model_abs[[mod]]))
  coeffs_cropland_abs <- summary(model_abs[[mod]])$coefficients %>%
    as.data.frame() %>%
    rownames_to_column("metric") %>%
    tibble() %>%
    mutate(gene = mod, .before = metric) %>%
    mutate(singular = sing) %>%
    rbind(coeffs_cropland_abs)
  
  print("")
  print("----------------------------------------------")
  print("----------------------------------------------")
  print("")
}

slopes <- list()
slopes$cropland$tpm <- coeffs_cropland_tpm_log %>%
  filter(metric == "cropland_fraction_2k_radius") %>%
  mutate(raw_significant = ifelse(`Pr(>|t|)` <= 0.05, "*", "n.s."),
         p_adjusted = p.adjust(`Pr(>|t|)`, method = "BH"),
         adjusted_significant = ifelse(p_adjusted <= 0.05, "*", "n.s."))
slopes$cropland$abs <- coeffs_cropland_abs %>%
  filter(metric == "cropland_fraction_2k_radius") %>%
  mutate(raw_significant = ifelse(`Pr(>|t|)` <= 0.05, "*", "n.s."),
         p_adjusted = p.adjust(`Pr(>|t|)`, method = "BH"),
         adjusted_significant = ifelse(p_adjusted <= 0.05, "*", "n.s."))

tpm_genes_of_particular_interest <- c("phosphoadenosine phosphosulfate reductase", "levanase", "PnuC-like nicotinamide mononucleotide transport")
abs_genes_of_particular_interest <- c("phosphoadenosine phosphosulfate reductase", "levanase")

summary(model_abs$`phosphoadenosine phosphosulfate reductase`)
summary(model_abs$levanase)

summary(model_tpm_log$`phosphoadenosine phosphosulfate reductase`)
summary(model_tpm_log$levanase)
summary(model_tpm_log$`PnuC-like nicotinamide mononucleotide transport`)

# Diagnostics (do this manually)
for (gopi in tpm_genes_of_particular_interest) {
  
  plot(model_tpm_log[[gopi]], which = 1)
  qqnorm(resid(model_tpm_log[[gopi]]))
  qqline(resid(model_tpm_log[[gopi]]))

  ranef_model <- ranef(model_tpm_log[[gopi]])$Hive_ID[[1]]
  hist(ranef_model, main = "Random Intercepts for Hive_ID", xlab = "Intercept")
  qqnorm(ranef_model)
  qqline(ranef_model)

  simulationOutput <- simulateResiduals(fittedModel = model_tpm_log[[gopi]])
  plot(simulationOutput)
  testUniformity(simulationOutput)
  testDispersion(simulationOutput)

}
for (gopi in abs_genes_of_particular_interest) {
  
  plot(model_abs[[gopi]], which = 1)
  qqnorm(resid(model_abs[[gopi]]))
  qqline(resid(model_abs[[gopi]]))
  
  ranef_model <- ranef(model_abs[[gopi]])$Hive_ID[[1]]
  hist(ranef_model, main = "Random Intercepts for Hive_ID", xlab = "Intercept")
  qqnorm(ranef_model)
  qqline(ranef_model)
  
  simulationOutput <- simulateResiduals(fittedModel = model_abs[[gopi]])
  plot(simulationOutput)
  testUniformity(simulationOutput)
  testDispersion(simulationOutput)
}
# 
# abs_raw_aggregated_cor_plot <- list()
# log_tpm_raw_aggregated_cor_plot <- list()
# for (gopi in genes_of_particular_interest) {
#   abs_raw_aggregated_cor_plot[[gopi]] <- test_tibble_abs[[gopi]] %>%
#     group_by(Country) %>%
#     mutate(mean_log_load = mean(log_load)) %>%
#     select(Country, cropland_fraction_2k_radius, mean_log_load) %>%
#     distinct() %>%
#     ggplot(aes(x = cropland_fraction_2k_radius, y = mean_log_load)) +
#     geom_point() +
#     geom_smooth(method = "glm", formula = y ~ x) +
#     stat_cor(method = "spearman") +
#     ggtitle(paste0(gopi, " - abs"))
#   
#   log_tpm_raw_aggregated_cor_plot[[gopi]] <- test_tibble_tpm_log[[gopi]] %>%
#     group_by(Country) %>%
#     mutate(mean_log_tpm = mean(log_tpm)) %>%
#     select(Country, cropland_fraction_2k_radius, mean_log_tpm) %>%
#     distinct() %>%
#     ggplot(aes(x = cropland_fraction_2k_radius, y = mean_log_tpm)) +
#     geom_point() +
#     geom_smooth(method = "glm", formula = y ~ x) +
#     stat_cor(method = "spearman") +
#     ggtitle(paste0(gopi, " - log_tpm"))
# }


# sulf gene and levanase never co-occur! 
# sulf gene and PnuC-like gene also never
# PnuC-like gene and levanase only once
grene_presence_on_contigs %>%
  select(contig, `phosphoadenosine phosphosulfate reductase`, levanase, `PnuC-like nicotinamide mononucleotide transport`) %>%
  mutate(
    occurrance_any = ifelse(`phosphoadenosine phosphosulfate reductase` == 1 |
                              levanase == 1 |
                              `PnuC-like nicotinamide mononucleotide transport` == 1,
                            1, 0),
    occurrance_sulf_and_levanase = ifelse(`phosphoadenosine phosphosulfate reductase` |
                                            levanase == 1,
                                          1, 0),
    occurrence_sulf_and_PnuC = ifelse(`phosphoadenosine phosphosulfate reductase` == 1 |
                                        `PnuC-like nicotinamide mononucleotide transport` == 1,
                                      1, 0),
    occurrence_levanase_and_PnuC = ifelse(levanase == 1 |
                                            `PnuC-like nicotinamide mononucleotide transport` == 1,
                                          1, 0)
    )

grene_presence_on_contigs %>%
  select(contig, all_of(tpm_genes_of_particular_interest)) %>%
  filter(`phosphoadenosine phosphosulfate reductase` == 1 | levanase == 1) %>%
  mutate(gopi = ifelse(levanase == 1, "levanase", "phosphoadenosine phosphosulfate reductase")) %>%
  select(contig, gopi) %>%
  left_join(., classification, by = "contig")



#### PESTICIDES:
#### ABS
pest_test_tibble_abs <- list()
pest_model_abs <- list()
coeffs_pesticides_abs <- tibble()
for (gopi in abs_genes_of_particular_interest) {
  for (item in c("Pesticides (total)", "Insecticides", "Herbicides", "Fungicides and Bactericides", "Plant Growth Regulators")) {
    pest_test_tibble_abs[[gopi]] <- FAOSTAT_added_data %>% 
      filter(Item == item) %>% 
      select(Country, est_use_in_2k_radius) %>%
      left_join(., test_tibble_abs[[gopi]], by = "Country")
    
    pest_model_abs[[gopi]][[item]] <- lmer(log_load ~ est_use_in_2k_radius + Season +
                               (1 | Hive_ID), data = pest_test_tibble_abs[[gopi]])
    coeffs_pesticides_abs <- summary(pest_model_abs[[gopi]][[item]])$coefficients %>%
      as.data.frame() %>%
      rownames_to_column("metric") %>%
      tibble() %>%
      mutate(gene = gopi, 
             Item = item,
             .before = metric) %>%
      mutate(singular = ifelse(isSingular(pest_model_abs[[gopi]][[item]]), TRUE, FALSE)) %>%
      rbind(coeffs_pesticides_abs)
  }
}

slopes$pesticides$abs <- coeffs_pesticides_abs %>%
  filter(metric == "est_use_in_2k_radius") %>%
  mutate(raw_significant = ifelse(`Pr(>|t|)` <= 0.05, "*", "n.s."),
         p_adjusted = p.adjust(`Pr(>|t|)`, method = "BH"),
         adjusted_significant = ifelse(p_adjusted <= 0.05, "*", "n.s."))
coeffs_pesticides_abs %>%
  filter(metric == "est_use_in_2k_radius",
         gene == "levanase") %>%
  mutate(raw_significant = ifelse(`Pr(>|t|)` <= 0.05, "*", "n.s."),
         p_adjusted = p.adjust(`Pr(>|t|)`, method = "BH"),
         adjusted_significant = ifelse(p_adjusted <= 0.05, "*", "n.s."))

summary(pest_model_abs$levanase$Herbicides)
summary(pest_model_abs$levanase$`Plant Growth Regulators`)

# Diagnostics
# mod <- "Herbicides"
# mod <- "Plant Growth Regulators"
# plot(pest_model_abs$levanase[[mod]], which = 1)
# qqnorm(resid(pest_model_abs$levanase[[mod]]))
# qqline(resid(pest_model_abs$levanase[[mod]]))
# simulationOutput <- simulateResiduals(fittedModel = pest_model_abs$levanase[[mod]])
# plot(simulationOutput)
# testUniformity(simulationOutput)
# testDispersion(simulationOutput)


herbs <- FAOSTAT_added_data %>%
  filter(str_detect(Item, "Herbicides ")) %>%
  distinct(Item) %>%
  unlist(use.names = FALSE)

herbs_test_tibble_abs <- list()
herbs_model_abs <- list()
coeffs_herbs_abs <- tibble()
for (item in herbs) {
  herbs_test_tibble_abs$levanase[[item]] <- FAOSTAT_added_data %>% 
    filter(Item == item) %>%
    select(Country, est_use_in_2k_radius) %>%
    left_join(., test_tibble_abs$levanase, by = "Country")

  herbs_model_abs$levanase[[item]] <- lmer(log_load ~ est_use_in_2k_radius + Season +
                                           (1 | Hive_ID), data = herbs_test_tibble_abs$levanase[[item]])
  
  coeffs_herbs_abs <- summary(herbs_model_abs$levanase[[item]])$coefficients %>%
    as.data.frame() %>%
    rownames_to_column("metric") %>%
    tibble() %>%
    mutate(gene = "levanase", 
           Item = item,
           .before = metric) %>%
    mutate(singular = ifelse(isSingular(herbs_model_abs$levanase[[item]]), TRUE, FALSE)) %>%
    rbind(coeffs_herbs_abs)
}

# Diagnostics
# mod <- "Herbicides – Phenoxy hormone products"
# mod <- "Herbicides – Triazines"
# mod <- "Herbicides – Amides"
# mod <- "Herbicides – Carbamates"
# mod <- "Herbicides – Dinitroanilines"
# mod <- "Herbicides – Urea derivates"
# mod <- "Herbicides – Sulfonyl ureas"
# mod <- "Herbicides – Bipiridils"
# mod <- "Herbicides – Uracil"
# mod <- "Herbicides - nes"
# plot(herbs_model_abs$levanase[[mod]], which = 1)
# qqnorm(resid(herbs_model_abs$levanase[[mod]]))
# qqline(resid(herbs_model_abs$levanase[[mod]]))
# simulationOutput <- simulateResiduals(fittedModel = herbs_model_abs$levanase[[mod]])
# plot(simulationOutput)
# testUniformity(simulationOutput)
# testDispersion(simulationOutput)

slopes$herbs$abs <- coeffs_herbs_abs %>%
  filter(metric == "est_use_in_2k_radius") %>%
  mutate(raw_significant = ifelse(`Pr(>|t|)` <= 0.05, "*", "n.s."),
         p_adjusted = p.adjust(`Pr(>|t|)`, method = "BH"),
         adjusted_significant = ifelse(p_adjusted <= 0.05, "*", "n.s."))

#### LOG TPM

pest_test_tibble_log_tpm <- list()
pest_model_log_tpm <- list()
coeffs_pesticides_log_tpm <- tibble()
for (gopi in tpm_genes_of_particular_interest) {
  for (item in c("Pesticides (total)", "Insecticides", "Herbicides", "Fungicides and Bactericides", "Plant Growth Regulators")) {
    pest_test_tibble_log_tpm[[gopi]] <- FAOSTAT_added_data %>% 
      filter(Item == item) %>% 
      select(Country, est_use_in_2k_radius) %>%
      left_join(., test_tibble_tpm_log[[gopi]], by = "Country")
    
    pest_model_log_tpm[[gopi]][[item]] <- lmer(log_tpm ~ est_use_in_2k_radius + Gut_part + Season +
                                             (1 | Hive_ID), data = pest_test_tibble_log_tpm[[gopi]])
    
    coeffs_pesticides_log_tpm <- summary(pest_model_log_tpm[[gopi]][[item]])$coefficients %>%
      as.data.frame() %>%
      rownames_to_column("metric") %>%
      tibble() %>%
      mutate(gene = gopi, 
             Item = item,
             .before = metric) %>%
      mutate(singular = ifelse(isSingular(pest_model_log_tpm[[gopi]][[item]]), TRUE, FALSE)) %>%
      rbind(coeffs_pesticides_log_tpm)
  }
}

# Diagnostics
# gene <- "phosphoadenosine phosphosulfate reductase"
# gene <- "levanase"
# gene <- "PnuC-like nicotinamide mononucleotide transport"
# mod <- "Pesticides (total)"
# mod <- "Insecticides"
# mod <- "Herbicides"
# mod <- "Fungicides and Bactericides"
# mod <- "Plant Growth Regulators"
# plot(pest_model_log_tpm$levanase[[mod]], which = 1)
# qqnorm(resid(pest_model_log_tpm$levanase[[mod]]))
# qqline(resid(pest_model_log_tpm$levanase[[mod]]))
# simulationOutput <- simulateResiduals(fittedModel = pest_model_log_tpm$levanase[[mod]])
# plot(simulationOutput)
# testUniformity(simulationOutput)
# testDispersion(simulationOutput)

slopes$pesticides$tpm <- coeffs_pesticides_log_tpm %>%
  filter(metric == "est_use_in_2k_radius") %>%
  mutate(raw_significant = ifelse(`Pr(>|t|)` <= 0.05, "*", "n.s."),
         p_adjusted = p.adjust(`Pr(>|t|)`, method = "BH"),
         adjusted_significant = ifelse(p_adjusted <= 0.05, "*", "n.s."))

herbs_and_fungs <- FAOSTAT_added_data %>%
  filter(str_detect(Item, "Herbicides ") |
           str_detect(Item, "Fung ")) %>%
  distinct(Item) %>%
  unlist(use.names = FALSE)

herbfung_tibble_log_tpm <- list()
herbfung_model_log_tpm <-list()
coeffs_herbs_and_fungs_log_tpm <- tibble()
for (item in herbs_and_fungs) {
  sulf <- "phosphoadenosine phosphosulfate reductase"
  herbfung_tibble_log_tpm[[sulf]][[item]] <- FAOSTAT_added_data %>% 
    filter(Item == item) %>% 
    select(Country, est_use_in_2k_radius) %>%
    left_join(., test_tibble_tpm_log[[sulf]], by = "Country")
  
  herbfung_model_log_tpm[[sulf]][[item]] <- lmer(log_tpm ~ est_use_in_2k_radius + Gut_part + Season +
                                            (1 | Hive_ID), data = herbfung_tibble_log_tpm[[sulf]][[item]])
  
  coeffs_herbs_and_fungs_log_tpm <- summary(herbfung_model_log_tpm[[sulf]][[item]])$coefficients %>%
    as.data.frame() %>%
    rownames_to_column("metric") %>%
    tibble() %>%
    mutate(gene = "phosphoadenosine phosphosulfate reductase", 
           Item = item,
           .before = metric) %>%
    mutate(singular = ifelse(isSingular(herbfung_model_log_tpm[[sulf]][[item]]), TRUE, FALSE)) %>%
    rbind(coeffs_herbs_and_fungs_log_tpm)
}

slopes$herbs_and_fungs$tpm <- coeffs_herbs_and_fungs_log_tpm %>%
  filter(metric == "est_use_in_2k_radius") %>%
  mutate(raw_significant = ifelse(`Pr(>|t|)` <= 0.05, "*", "n.s."),
         p_adjusted = p.adjust(`Pr(>|t|)`, method = "BH"),
         adjusted_significant = ifelse(p_adjusted <= 0.05, "*", "n.s."))


## THE STORY:

# Checking the 15 metabolism genes for correlation with cropland in area, 
# 12 of which had valid input data (loggable)
slopes$cropland$tpm
# Sulf gene is significantly (negatively) correlated 

# Checking if sulf gene correlates with pesticide use in the surrounding cropland
slopes$pesticides$tpm %>% 
  filter(gene == "phosphoadenosine phosphosulfate reductase",
         Item == "Pesticides (total)")
# It does.

# Checking of sulf gene correlates with any of the larger pesticide groups
slopes$pesticides$tpm %>% 
  filter(gene == "phosphoadenosine phosphosulfate reductase",
         Item != "Pesticides (total)") %>%
  mutate(p_adjusted = p.adjust(`Pr(>|t|)`))
# Total pesticide use, Herbicide use, Fungicide use and plant growth regulator use correlate significantly (insecticides don't).

# Going deeper into the different herbicides and fungicides (there are no sub-categories of plant growth regs)
slopes$herbs_and_fungs$tpm
# A lot of them are significant!


