# library(stageR)
library(lme4)
library(lmerTest)
library(DHARMa)
library(ggrepel)
library(tidyverse)
library(patchwork)
# library(car)
# library(emmeans)
# library(ggpubr)
# library(gMCPLite)
library(gMCP)

forest_plot <- function(tbl, axis_name = NULL, plot_title = NULL) {
  tbl %>%
    ggplot(aes(x = Estimate, y = axis_labels)) +
    geom_pointrange(aes(xmin = Estimate - `Std. Error`,
                        xmax = Estimate + `Std. Error`,
                        alpha = p_adjusted <= 0.05),  # Adjust alpha for significance
                    color = "black",
                    size = 0.7) +
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.3)) +
    # theme_minimal() +
    
    # Place p_adjusted labels at the far right with color mapped to significance
    # geom_text(aes(x = Inf, 
    #               label = sprintf("%.3f", p_adjusted), 
    #               color = p_adjusted <= 0.05), 
    #           hjust = -0.1,    # adjust horizontal justification as needed
    #           size = 3) +
    # scale_color_manual(values = c("TRUE" = "black", "FALSE" = "gray")) +
    geom_text_repel(aes(label = sprintf("%.3f", p_adjusted),
                        color = p_adjusted <= 0.05),
                    direction = "y",
                    nudge_y = 0.1,
                    segment.size = 0,
                    box.padding = 0.5) +
    scale_color_manual(values = c("TRUE" = "black", "FALSE" = "gray")) +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(y = axis_name,
         x = "Estimated slope (log)") +
    ggtitle(plot_title)
}


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
  left_join(cropland_fraction[c("Country", "cropland_fraction_2k_radius")], ., by = "Country") %>%
  select_if(~ !any(is.na(.)))

cropland_and_FAO_ranks <- cropland_and_FAO %>%
  mutate(across(-Country, ~ rank(.x, ties.method = "average")))

# cropland_and_FAO %>% 
#   ggplot(aes(x = cropland_fraction_2k_radius, y = `Plant Growth Regulators`)) +
#   geom_point() +
#   geom_smooth(method = "glm", formula = y ~ x) +
#   stat_cor(method = "spearman")

metadata <- readRDS("output/R/R_variables/metadata.RDS") %>%
  mutate(Hive_ID = as.character(Hive_ID))
classification <- readRDS("output/R/R_variables/classification.RDS")

phage_tpm <- read.csv("output/R/relative_abundance/phage_tpm.csv") %>%
  tibble()

# absolute_counts <- read.csv("output/R/absolute_counts.csv") %>%
#   tibble()

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
  select(-contig) %>%
  distinct()

# gene_load <- grene_presence_on_contigs %>%
#   pivot_longer(-contig, names_to = "gene", values_to = "present") %>%
#   inner_join(., absolute_counts, by = "contig") %>%
#   mutate(across(-c(contig, gene, present), ~ .x * present)) %>%
#   select(-present) %>%
#   pivot_longer(-c(contig, gene), names_to = "Sample_ID", values_to = "load") %>%
#   group_by(gene, Sample_ID) %>%
#   mutate(load = sum(load)) %>%
#   ungroup() %>%
#   select(-contig) %>%
#   distinct()


genes_of_interest <- unique(gene_tpm_long$Sample_ID$gene)

gene_tpm_long$Sample_ID

test_tibble_tpm_log <- list()
model_tpm_log <- list()

model_anova <- list()
pw_results <- list()
emm_plots <- list()
emm_df <- list()

# goi <- "phosphoadenosine phosphosulfate reductase" # sing # CALCULATED
# goi <- "levanase" # CALCULATED
# goi <- "ribosomal protein S6 glutaminyl transferase"
# goi <- "glucosyltransferase" # CALCULATED
# goi <- "glutamine amidotransferase"
# goi <- "porphyrin biosynthesis" # Warining
# goi <- "aerobic cobaltochelatase CobT subunit" # Warning
# goi <- "chitinase" # CALCULATED
# goi <- "dTDP-4-dehydrorhamnose 3"
# goi <- "nicotinamide-nucleotide adenylyltransferase" # Sing
# goi <- "PnuC-like nicotinamide mononucleotide transport" # CALCULATED
# goi <- "NrdD-like anaerobic ribonucleotide reductase large subunit" #Column drop, sing
# goi <- "QueD-like  6-pyruvoyl-tetrahydropterin synthase" # No model
# goi <- "QueC-like queuosine biosynthesis" # No model
# goi <- "QueE-like  radical SAM domain" # No model

for (goi in genes_of_interest) {
  
  test_tibble_tpm_log[[goi]] <- gene_tpm_long$Sample_ID %>%
    filter(gene == goi) %>%
    left_join(., metadata[c("Sample_ID", "Country", "Hive_ID", "Season", "Gut_part")], by = "Sample_ID") %>%
    distinct() %>%
    # left_join(., cropland_fraction[c("Country", "cropland_fraction_2k_radius")], by = "Country") %>%
    left_join(., cropland_and_FAO, by = "Country") %>%
    # left_join(., cropland_and_FAO_ranks, by = "Country") %>%
    mutate(Gut_part = factor(Gut_part, levels = c("rec", "ile", "mid"))) %>%
    mutate(log_tpm = log10(tpm)) %>%
    filter(!is.infinite(log_tpm))
  
  countries_left <- test_tibble_tpm_log[[goi]] %>%
    distinct(Country) %>%
    nrow()
  
  if (countries_left == 8 ) {
    model_tpm_log[[goi]] <- lmer(log_tpm ~ cropland_fraction_2k_radius + Gut_part + Season +
                                   (1 | Hive_ID), data = test_tibble_tpm_log[[goi]])
    
    # model_tpm_log[[goi]] <- lmer(log_tpm ~ cropland_fraction_2k_radius * Gut_part +
    #        cropland_fraction_2k_radius * Season + (1 | Hive_ID),
    #      data = test_tibble_tpm_log[[goi]])
    
    
    # As categorical
    # model_tpm_log[[goi]] <- lmer(log_tpm ~ Country + Gut_part + Season +
    #                                (1 | Hive_ID), data = test_tibble_tpm_log[[goi]])
    # 
    # summary(model_tpm_log[[goi]])

    # 
    # model_anova[[goi]] <- Anova(model_tpm_log[[goi]], type = "III")
    # # model_full <- lmer(log_tpm ~ Country + Gut_part + Season + (1 | Hive_ID),
    # #                    data = test_tibble_tpm_log[[goi]])
    # # model_reduced <- lmer(log_tpm ~ Gut_part + Season + (1 | Hive_ID),
    # #                       data = test_tibble_tpm_log[[goi]])
    # # 
    # # 
    # # anova(model_reduced, model_full)
    # 
    # emm_country <- emmeans(model_tpm_log[[goi]], "Country")
    # pw_results[[goi]] <- pairs(emm_country)
    # 
    # emm_df[[goi]] <- as.data.frame(emm_country)
    # emm_plots[[goi]] <- ggplot(emm_df[[goi]], aes(x = Country, y = emmean)) +
    #   geom_point() +
    #   geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2) +
    #   labs(x = "Country", y = "Estimated mean log(tpm)",
    #        title = paste0("Estimated marginal means by Country - ", goi)) +
    #   theme_minimal()
    
  } else {
    # print(paste0(goi, " - model for log tpm not possible. Hive_IDs left: ", hive_ids, ". Samples left: ", nrow(test_tibble_tpm_log[[goi]]), "."))
    print(paste0(goi, " - no model calculated. Only ", countries_left, " countries left."))
  }
}


significant_in_anova <- c("phosphoadenosine phosphosulfate reductase", "levanase")

items_of_interest <- c("cropland_fraction_2k_radius", "Pesticides (total)", "Insecticides",
                       "Herbicides", "Fungicides and Bactericides")
cropland_and_FAO %>% colnames()

# 
# trend_tibble <- list()
# trend_plot <- list()
# trend_plot_ranked <- list()
# for (goi in significant_in_anova) {
#   # for (item in cropland_and_FAO %>% select(-Country) %>% colnames()) {
#   for (item in items_of_interest) {
#     trend_tibble[[goi]][[item]] <- emm_df[[goi]] %>%
#       as_tibble() %>%
#       full_join(., cropland_and_FAO, by = "Country") %>%
#       select(Country, emmean, all_of(item)) %>%
#       pivot_longer(-Country) %>%
#       group_by(name) %>%
#       mutate(max_abs = max(abs(value))) %>%
#       ungroup() %>%
#       mutate(scaling_factor = max(max_abs) / min(max_abs)) %>%
#       mutate(value_scaled = ifelse(name == item, value * scaling_factor, -value)) %>%
#       mutate(abs_value = abs(value)) %>%
#       group_by(name) %>%
#       mutate(rank_value = rank(abs_value))
#         
#     scaling_factor <- unique(trend_tibble[[goi]][[item]]$scaling_factor)
#       
#     trend_plot[[goi]][[item]] <- ggplot(trend_tibble[[goi]][[item]], aes(x = Country, y = value_scaled, fill = name)) +
#       geom_col(position = "dodge") +
#       ggtitle(paste0(goi, " - ", item)) +
#       scale_y_continuous(
#         name = "emmean",
#         sec.axis = sec_axis(~ . / scaling_factor , name = item))
#     
#     trend_plot_ranked[[goi]][[item]] <- ggplot(trend_tibble[[goi]][[item]], aes(x = Country, y = rank_value, fill = name)) +
#       geom_col(position = "dodge") +
#       ggtitle(paste0(goi, " - ", item))
#     
#   }
# }

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

slopes <- list()
slopes$cropland$tpm <- coeffs_cropland_tpm_log %>%
  filter(metric == "cropland_fraction_2k_radius") %>%
  mutate(raw_significant = ifelse(`Pr(>|t|)` <= 0.05, "*", "n.s."),
         p_adjusted = p.adjust(`Pr(>|t|)`, method = "BH"),
         adjusted_significant = ifelse(p_adjusted <= 0.05, "*", "n.s."))

tpm_genes_of_particular_interest <- c("phosphoadenosine phosphosulfate reductase", "levanase", "PnuC-like nicotinamide mononucleotide transport")

summary(model_tpm_log$`phosphoadenosine phosphosulfate reductase`)
summary(model_tpm_log$levanase)
summary(model_tpm_log$`PnuC-like nicotinamide mononucleotide transport`)

# Diagnostics (do this manually)

# gopi <- "phosphoadenosine phosphosulfate reductase"
# gopi <- "levanase"
# gopi <- "PnuC-like nicotinamide mononucleotide transport"
# for (gopi in tpm_genes_of_particular_interest) {
#   
#   plot(model_tpm_log[[gopi]], which = 1)
#   qqnorm(resid(model_tpm_log[[gopi]]))
#   qqline(resid(model_tpm_log[[gopi]]))
#   
#   ranef_model <- ranef(model_tpm_log[[gopi]])$Hive_ID[[1]]
#   hist(ranef_model, main = "Random Intercepts for Hive_ID", xlab = "Intercept")
#   qqnorm(ranef_model)
#   qqline(ranef_model)
#   
#   simulationOutput <- simulateResiduals(fittedModel = model_tpm_log[[gopi]])
#   plot(simulationOutput)
#   testUniformity(simulationOutput)
#   testDispersion(simulationOutput)
#   
# }

# grene_presence_on_contigs %>%
#   select(contig, `phosphoadenosine phosphosulfate reductase`, levanase, `PnuC-like nicotinamide mononucleotide transport`) %>%
#   mutate(
#     occurrance_any = ifelse(`phosphoadenosine phosphosulfate reductase` == 1 |
#                               levanase == 1 |
#                               `PnuC-like nicotinamide mononucleotide transport` == 1,
#                             1, 0),
#     occurrance_sulf_and_levanase = ifelse(`phosphoadenosine phosphosulfate reductase` |
#                                             levanase == 1,
#                                           1, 0),
#     occurrence_sulf_and_PnuC = ifelse(`phosphoadenosine phosphosulfate reductase` == 1 |
#                                         `PnuC-like nicotinamide mononucleotide transport` == 1,
#                                       1, 0),
#     occurrence_levanase_and_PnuC = ifelse(levanase == 1 |
#                                             `PnuC-like nicotinamide mononucleotide transport` == 1,
#                                           1, 0)
#   )
# 
# grene_presence_on_contigs %>%
#   select(contig, all_of(tpm_genes_of_particular_interest)) %>%
#   filter(`phosphoadenosine phosphosulfate reductase` == 1 | levanase == 1) %>%
#   mutate(gopi = ifelse(levanase == 1, "levanase", "phosphoadenosine phosphosulfate reductase")) %>%
#   select(contig, gopi) %>%
#   left_join(., classification, by = "contig")
# 


#### PESTICIDES:
#### LOG TPM
total_pest_test_tibble_log_tpm <- list()
total_pest_model_log_tpm <- list()
coeffs_total_pest_log_tpm <- tibble()
for (gopi in tpm_genes_of_particular_interest) {
  for (item in "Pesticides (total)") {
    # total_pest_test_tibble_log_tpm[[gopi]] <- FAOSTAT_added_data %>% 
    #   filter(Item == item) %>% 
    #   select(Country, est_use_in_2k_radius) %>%
    #   left_join(., test_tibble_tpm_log[[gopi]], by = "Country")
    
    # total_pest_model_log_tpm[[gopi]][[item]] <- lmer(log_tpm ~ est_use_in_2k_radius + Gut_part + Season +
    #                                              (1 | Hive_ID), data = total_pest_test_tibble_log_tpm[[gopi]])
    
    temp_test_tibble <- test_tibble_tpm_log[[gopi]] %>% 
      rename(est_use_in_2k_radius = all_of(item))
    total_pest_model_log_tpm[[gopi]][[item]] <- lmer(log_tpm ~ est_use_in_2k_radius + Gut_part + Season +
                                                       (1 | Hive_ID), data = temp_test_tibble)
    
    coeffs_total_pest_log_tpm <- summary(total_pest_model_log_tpm[[gopi]][[item]])$coefficients %>%
      as.data.frame() %>%
      rownames_to_column("metric") %>%
      tibble() %>%
      mutate(gene = gopi, 
             Item = item,
             .before = metric) %>%
      mutate(singular = ifelse(isSingular(total_pest_model_log_tpm[[gopi]][[item]]), TRUE, FALSE)) %>%
      rbind(coeffs_total_pest_log_tpm)
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
# plot(total_pest_model_log_tpm$levanase[[mod]], which = 1)
# qqnorm(resid(total_pest_model_log_tpm$levanase[[mod]]))
# qqline(resid(total_pest_model_log_tpm$levanase[[mod]]))
# simulationOutput <- simulateResiduals(fittedModel = total_pest_model_log_tpm$levanase[[mod]])
# plot(simulationOutput)
# testUniformity(simulationOutput)
# testDispersion(simulationOutput)

slopes$pesticides$tpm <- coeffs_total_pest_log_tpm %>%
  filter(metric == "est_use_in_2k_radius") %>%
  mutate(raw_significant = ifelse(`Pr(>|t|)` <= 0.05, "*", "n.s."),
         p_adjusted = p.adjust(`Pr(>|t|)`, method = "BH"),
         adjusted_significant = ifelse(p_adjusted <= 0.05, "*", "n.s."))


pest_groups_test_tibble_log_tpm <- list()
pest_groups_model_log_tpm <- list()
coeffs_pest_groups_log_tpm <- tibble()
for (gopi in tpm_genes_of_particular_interest) {
  for (item in c("Insecticides", "Herbicides", "Fungicides and Bactericides", "Plant Growth Regulators")) {
    # pest_groups_test_tibble_log_tpm[[gopi]] <- FAOSTAT_added_data %>% 
    #   filter(Item == item) %>% 
    #   select(Country, est_use_in_2k_radius) %>%
    #   left_join(., test_tibble_tpm_log[[gopi]], by = "Country")
    # 
    # pest_groups_model_log_tpm[[gopi]][[item]] <- lmer(log_tpm ~ est_use_in_2k_radius + Gut_part + Season +
    #                                              (1 | Hive_ID), data = pest_groups_test_tibble_log_tpm[[gopi]])
    
    temp_test_tibble <- test_tibble_tpm_log[[gopi]] %>% 
      rename(est_use_in_2k_radius = all_of(item))
    pest_groups_model_log_tpm[[gopi]][[item]] <- lmer(log_tpm ~ est_use_in_2k_radius + Gut_part + Season +
                                                       (1 | Hive_ID), data = temp_test_tibble)
    
    coeffs_pest_groups_log_tpm <- summary(pest_groups_model_log_tpm[[gopi]][[item]])$coefficients %>%
      as.data.frame() %>%
      rownames_to_column("metric") %>%
      tibble() %>%
      mutate(gene = gopi, 
             Item = item,
             .before = metric) %>%
      mutate(singular = ifelse(isSingular(pest_groups_model_log_tpm[[gopi]][[item]]), TRUE, FALSE)) %>%
      rbind(coeffs_pest_groups_log_tpm)
  }
}

slopes$pest_groups$tpm <- coeffs_pest_groups_log_tpm %>%
  filter(metric == "est_use_in_2k_radius") %>%
  mutate(raw_significant = ifelse(`Pr(>|t|)` <= 0.05, "*", "n.s."),
         p_adjusted = p.adjust(`Pr(>|t|)`, method = "BH"),
         adjusted_significant = ifelse(p_adjusted <= 0.05, "*", "n.s."))



herbs <- test_tibble_tpm_log[[gopi]] %>% 
  colnames() %>% 
  tibble() %>%
  filter(str_detect(., "Herbicides ")) %>%
  distinct(.) %>%
  unlist(use.names = FALSE)
fungs <- test_tibble_tpm_log[[gopi]] %>% 
  colnames() %>% 
  tibble() %>%
  filter(str_detect(., "Fung & Bact ")) %>%
  distinct(.) %>%
  unlist(use.names = FALSE)
insects <- test_tibble_tpm_log[[gopi]] %>% 
  colnames() %>% 
  tibble() %>%
  filter(str_detect(., "Insecticides ")) %>%
  distinct(.) %>%
  unlist(use.names = FALSE)

gene_pest_combos <- list("phosphoadenosine phosphosulfate reductase" = c(herbs, fungs),
     "levanase" = herbs,
     "PnuC-like nicotinamide mononucleotide transport" = c(herbs, fungs, insects)
     )


specific_pests_tibble_log_tpm <- list()
specific_pests_model_log_tpm <-list()
coeffs_specific_pests_log_tpm <- tibble()
for (gopi in tpm_genes_of_particular_interest) {
  for (item in gene_pest_combos[[gopi]]) {
    # specific_pests_tibble_log_tpm[[gopi]][[item]] <- FAOSTAT_added_data %>% 
    #   filter(Item == item) %>%
    #   select(Country, est_use_in_2k_radius) %>%
    #   left_join(., test_tibble_tpm_log[[gopi]], by = "Country")
    # 
    # specific_pests_model_log_tpm[[gopi]][[item]] <- lmer(log_tpm ~ est_use_in_2k_radius + Gut_part + Season +
    #        (1 | Hive_ID), data = specific_pests_tibble_log_tpm[[gopi]][[item]])
    
    temp_test_tibble <- test_tibble_tpm_log[[gopi]] %>% 
      rename(est_use_in_2k_radius = all_of(item))
    specific_pests_model_log_tpm[[gopi]][[item]] <- lmer(log_tpm ~ est_use_in_2k_radius + Gut_part + Season +
                                                        (1 | Hive_ID), data = temp_test_tibble)
    
    coeffs_specific_pests_log_tpm <- summary(specific_pests_model_log_tpm[[gopi]][[item]])$coefficients %>%
      as.data.frame() %>%
      rownames_to_column("metric") %>%
      tibble() %>%
      mutate(gene = gopi, 
             Item = item,
             .before = metric) %>%
      mutate(singular = ifelse(isSingular(specific_pests_model_log_tpm[[gopi]][[item]]), TRUE, FALSE)) %>%
      rbind(coeffs_specific_pests_log_tpm)
    
  }
}

slopes$specific_pests$tpm <- coeffs_specific_pests_log_tpm %>%
  filter(metric == "est_use_in_2k_radius") %>%
  mutate(raw_significant = ifelse(`Pr(>|t|)` <= 0.05, "*", "n.s."),
         p_adjusted = p.adjust(`Pr(>|t|)`, method = "BH"),
         adjusted_significant = ifelse(p_adjusted <= 0.05, "*", "n.s."))


#########
# 
# herbs_and_fungs <- FAOSTAT_added_data %>%
#   filter(str_detect(Item, "Herbicides ") |
#            str_detect(Item, "Fung ")) %>%
#   distinct(Item) %>%
#   unlist(use.names = FALSE)
# 
# herbfung_tibble_log_tpm <- list()
# herbfung_model_log_tpm <-list()
# coeffs_herbs_and_fungs_log_tpm <- tibble()
# # "Herbicides – Uracil" predictor on different scales?
# for (item in herbs_and_fungs) {
#   sulf <- "phosphoadenosine phosphosulfate reductase"
#   herbfung_tibble_log_tpm[[sulf]][[item]] <- FAOSTAT_added_data %>% 
#     filter(Item == item) %>% 
#     select(Country, est_use_in_2k_radius) %>%
#     left_join(., test_tibble_tpm_log[[sulf]], by = "Country")
#   
#   herbfung_model_log_tpm[[sulf]][[item]] <- lmer(log_tpm ~ est_use_in_2k_radius + Gut_part + Season +
#                                                    (1 | Hive_ID), data = herbfung_tibble_log_tpm[[sulf]][[item]])
#   
#   coeffs_herbs_and_fungs_log_tpm <- summary(herbfung_model_log_tpm[[sulf]][[item]])$coefficients %>%
#     as.data.frame() %>%
#     rownames_to_column("metric") %>%
#     tibble() %>%
#     mutate(gene = "phosphoadenosine phosphosulfate reductase", 
#            Item = item,
#            .before = metric) %>%
#     mutate(singular = ifelse(isSingular(herbfung_model_log_tpm[[sulf]][[item]]), TRUE, FALSE)) %>%
#     rbind(coeffs_herbs_and_fungs_log_tpm)
# }
# 
# slopes$herbs_and_fungs$tpm <- coeffs_herbs_and_fungs_log_tpm %>%
#   filter(metric == "est_use_in_2k_radius") %>%
#   mutate(raw_significant = ifelse(`Pr(>|t|)` <= 0.05, "*", "n.s."),
#          p_adjusted = p.adjust(`Pr(>|t|)`, method = "BH"),
#          adjusted_significant = ifelse(p_adjusted <= 0.05, "*", "n.s."))
# 

## THE STORY:

layered_p_values <- list()
# layered_p_values$cropland <- slopes$cropland$tpm$p_adjusted %>%
layered_p_values$cropland <- slopes$cropland$tpm$`Pr(>|t|)` %>%
  setNames(slopes$cropland$tpm$gene)
layered_p_values$pesticides <- slopes$pesticides$tpm %>%
  filter(gene == "phosphoadenosine phosphosulfate reductase") %>%
  select(Item, `Pr(>|t|)`) %>%
  deframe()
layered_p_values$pest_groups <- slopes$pest_groups$tpm %>%
  filter(gene == "phosphoadenosine phosphosulfate reductase") %>%
  select(Item, `Pr(>|t|)`) %>%
  deframe()
layered_p_values$specific_pests <- slopes$specific_pests$tpm %>%
  filter(gene == "phosphoadenosine phosphosulfate reductase") %>%
  select(Item, `Pr(>|t|)`) %>%
  deframe()

layered_p_adjustments <- function(p_value_list = layered_p_values) {

  hypotheses <- c(names(p_value_list$cropland),
                  names(p_value_list$pesticides),
                  names(p_value_list$pest_groups),
                  names(p_value_list$specific_pests)
                  )

  herbs_in_hypoths <- names(p_value_list$specific_pests) %>%
    tibble() %>%
    filter(str_detect(., "Herbicides ")) %>%
    distinct(.) %>%
    unlist(use.names = FALSE)
  fungs_in_hypoths <- names(p_value_list$specific_pests) %>%
    tibble() %>%
    filter(str_detect(., "Fung & Bact ")) %>%
    distinct(.) %>%
    unlist(use.names = FALSE)

  # Define the initial weights.
  # Here we allocate alpha equally to the Stage 1 tests (genes) and zero to the rest.
  # When the hypothesis of a gene-test is rejected, its weight will be transferred.
  w <- rep(0, length(hypotheses))
  names(w) <- hypotheses
  w[names(p_value_list$cropland)] <- 1 / length(p_value_list$cropland)

  # Create an empty transition matrix.
  transitions <- matrix(0, nrow = length(hypotheses), ncol = length(hypotheses),
                        dimnames = list(hypotheses, hypotheses))
  # Define transitions from Stage 1 (genes) to Stage 2 (total pesticide test).
  # Each gene for which the hypothesis is rejected passes its entire weight to the Stage 2 test.
  for(g in names(p_value_list$cropland)){
    transitions[g, "Pesticides (total)"] <- 1
  }
  # Define transitions from Stage 2 (TotalPest) to Stage 3 (pesticide groups).
  # If the TotalPest test is rejected, split its weight equally among all 4 group tests.
  transitions["Pesticides (total)", names(p_value_list$pest_groups)] <- rep(1/length(p_value_list$pest_groups), length(p_value_list$pest_groups))
  # Define transitions from Stage 3 to Stage 4.
  # Only groups with further tests pass on weight.
  # For Group2 (assumed to have 6 specific pesticide tests):
  transitions["Fungicides and Bactericides", fungs_in_hypoths] <- rep(1/length(fungs_in_hypoths), length(fungs_in_hypoths))
  # For Group3 (assumed to have 7 specific pesticide tests):
  transitions["Herbicides", herbs_in_hypoths] <- rep(1/length(herbs_in_hypoths), length(herbs_in_hypoths))
  # Groups that are terminal (e.g., Group1 and Group4) need no transitions.
  
  # Read the transition matrix line by line, where the row name is the test in 
  # question and the column name is the test in the next layer of testing.
  # A 0 means that this test in the row name is not propagted to the test in the 
  # col name. A number different to 0 is the weight this test (row name) will 
  # have on the test in the next layer (col name). The weights in each line 
  # should sum up to 1. In an easy setting, equally divide the weight to the
  # test of the next layer.

  graph <- matrix2graph(m = transitions, weights = w)

  p_raw <- c(p_value_list$cropland, p_value_list$pesticides, p_value_list$pest_groups, p_value_list$specific_pests)
  # Optionally, inspect the graph.
  print(graph)
  
  cor_mat <- matrix(as.numeric(NA), nrow = length(hypotheses), ncol = length(hypotheses),
                                     dimnames = list(hypotheses, hypotheses))
  cor_mat[row(cor_mat)==col(cor_mat)] = 1

  result1 <- graphTest(pvalues = p_raw, graph = graph, cr = cor_mat)
  result2 <- gMCP(pvalues = p_raw, graph = graph, correlation = cor_mat)
  print(result)
  adjusted_p <- result@adjPValues

}


# Checking the 15 metabolism genes for correlation with cropland in area, 
# 12 of which had valid input data (loggable) and 5 had presence in all countries
slopes$cropland$tpm
# Sulf gene is significantly (negatively) correlated 

# Checking if sulf gene correlates with pesticide use in the surrounding cropland
focused_slopes_pesticides <- slopes$pesticides$tpm %>% 
  filter(gene == "phosphoadenosine phosphosulfate reductase")
focused_slopes_pesticides
# It does.

# Checking of sulf gene correlates with any of the larger pesticide groups
focused_slopes_pest_groups <- slopes$pest_groups$tpm %>% 
  filter(gene == "phosphoadenosine phosphosulfate reductase") %>%
  mutate(p_adjusted = p.adjust(`Pr(>|t|)`),
         adjusted_significant = ifelse(p_adjusted <= 0.05, "*", "n.s."))
focused_slopes_pest_groups
# Total pesticide use, Herbicide use, Fungicide use and plant growth regulator use correlate significantly (insecticides don't).

# Going deeper into the different herbicides and fungicides (there are no sub-categories of plant growth regs)
focused_slopes_specific_pests <- slopes$specific_pests$tpm %>%
  filter(gene == "phosphoadenosine phosphosulfate reductase")
focused_slopes_specific_pests
# A lot of them are significant!

focused_forests <- list()
focused_forests$cropland <- slopes$cropland$tpm %>%
  mutate(axis_labels = gene) %>%
  forest_plot(plot_title = "Cropland fraction in 2 km radius")
focused_forests$total_pest <- focused_slopes_pesticides %>%
  mutate(axis_labels = gene) %>%
  forest_plot(plot_title = "Total pesticide use in 2 km radius")
focused_forests$pest_groups <- focused_slopes_pest_groups %>%
  mutate(axis_labels = Item) %>%
  forest_plot(plot_title = "phosphoadenosine phosphosulfate reductase - Use of pesticide types in 2 km radius")
focused_forests$specific_pests <- focused_slopes_specific_pests %>% 
  filter(!Item %in% c("Herbicides – Uracil", "Herbicides – Sulfonyl ureas", "Herbicides – Bipiridils")) %>%
  mutate(axis_labels = Item) %>%
  forest_plot(plot_title = "phosphoadenosine phosphosulfate reductase -  Use of specific pesticides in 2 km radius")

focused_wrap_ratios <- c(
  length(focused_forests$cropland[["data"]][["axis_labels"]]),
  length(focused_forests$total_pest[["data"]][["axis_labels"]]),
  length(focused_forests$pest_groups[["data"]][["axis_labels"]]),
  length(focused_forests$specific_pests[["data"]][["axis_labels"]]
         ))
focused_forests_wrap <- wrap_plots(focused_forests, ncol = 1, 
                                   heights = focused_wrap_ratios, 
                                   axes = "collect")

## The side story
unfocused_forests <- list()
unfocused_forests$corpland <- focused_forests$cropland
unfocused_forests$total_pest <- slopes$pesticides$tpm %>%
  # filter(gene != "phosphoadenosine phosphosulfate reductase") %>%
  mutate(axis_labels = gene) %>%
  forest_plot(plot_title = "Total pesticide use in 2 km radius")
unfocused_forests$pest_groups <- slopes$pest_groups$tpm %>%
  # filter(gene != "phosphoadenosine phosphosulfate reductase") %>%
  mutate(axis_labels = paste0(gene, "; ", Item)) %>%
  forest_plot(plot_title = "Use of pesticide types in 2 km radius")
unfocused_forests$specific_pests <- slopes$specific_pests$tpm %>%
  # filter(gene != "phosphoadenosine phosphosulfate reductase") %>%
  mutate(axis_labels = paste0(gene, "; ", Item)) %>%
  forest_plot(plot_title = "Use of specific pesticides in 2 km radius")
unfocused_wrap_ratios <- c(
  length(unfocused_forests$corpland[["data"]][["axis_labels"]]),
  length(unfocused_forests$total_pest[["data"]][["axis_labels"]]),
  length(unfocused_forests$pest_groups[["data"]][["axis_labels"]]),
  length(unfocused_forests$specific_pests[["data"]][["axis_labels"]]
  ))
unfocused_forests_wrap <- wrap_plots(unfocused_forests, ncol = 1, 
                                   heights = unfocused_wrap_ratios, 
                                   axes = "collect")

# Save files
# system("mkdir -p output/R/gene_content/landuse/")
# ggsave("output/R/gene_content/landuse/focused_forests.pdf",
#        focused_forests_wrap, width = 15, height = 14)
# ggsave("output/R/gene_content/landuse/unfocused_forests.pdf",
#        unfocused_forests_wrap, width = 15, height = 35)
# write_delim(gene_tpm_long$Sample_ID, "output/R/gene_content/landuse/gene_tpm.tsv", 
#             delim = "\t")
