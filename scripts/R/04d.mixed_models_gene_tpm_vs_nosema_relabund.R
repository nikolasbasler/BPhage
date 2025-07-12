library(lme4)
library(lmerTest)
library(ggrepel)
library(forcats)
library(tidyverse)
library(patchwork)
library(tidytext)
library(ggtext)


source("scripts/R/helpers/mixed_helpers.R")

metadata <- readRDS("output/R/R_variables/metadata.RDS") %>%
  mutate(Hive_ID = as.character(Hive_ID))
classification <- readRDS("output/R/R_variables/classification.RDS")

nosema_relabund <- read.delim("output/nosema_mapped_counts_all.tsv") %>% 
  tibble() %>%
  mutate(nosema_relabund = mapped_to_Nosema / Hostout_R1_plus_R2) %>%
  select(Sample_ID, nosema_relabund)

phage_tpm <- read.csv("output/R/relative_abundance/phage_tpm.csv") %>%
  tibble()

phold_predictions_with_extensions <- read.delim("output/R/gene_content/phold_predictions_with_extensions_bphage_renamed_genes.tsv")

kegg_mapping <- read.delim("data/kegg_mapping.tsv", colClasses = "character") %>%
  tibble()

kegg_and_phold <- kegg_mapping %>%
  left_join(., phold_predictions_with_extensions[c("cds_id", "phrog", "function.", "product")], by = "cds_id")

CDSs_with_metabolism_kegg <- kegg_and_phold %>%
  filter(Pathway_category == "Metabolism" | 
           product %in% c("Chitinase", "GATase", "PnuC")) %>%
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

genes_of_interest <- c("Chitinase",
                       "Glucosyltransferase",
                       "Levanase",
                       "PAPS reductase", 
                       "PnuC")

coeffs <- list()

#####
# Gene TPM vs Nosema relabund
test_tibble_log_tpm <- list()
model_tpm <- list()
coeffs$tpm <- tibble()
for (goi in genes_of_interest) {
  
  test_tibble_log_tpm[[goi]] <- gene_tpm %>%
    filter(gene == goi) %>%
    left_join(., metadata[c("Sample_ID", "Country", "Hive_ID", "Season", "Gut_part")], by = "Sample_ID") %>%
    distinct() %>%
    left_join(., nosema_relabund, by = "Sample_ID") %>%
    mutate(Gut_part = factor(Gut_part, levels = c("rec", "ile", "mid"))) %>%
    mutate(log_tpm = log10(tpm),
           log_Vairimorpha_relabund = log10(nosema_relabund)) %>%
    filter(!is.infinite(log_tpm),
           !is.infinite(log_Vairimorpha_relabund))
  
  model_tpm[[goi]] <- lmer(log_tpm ~ log_Vairimorpha_relabund + Gut_part + Season +
                                             (1 | Hive_ID ), data = test_tibble_log_tpm[[goi]])
  
  has_convergence_issues <- FALSE
  messages <- model_tpm[[goi]]@optinfo$conv$lme4$messages
  if (is.null(messages)) {
    has_convergence_issues <- FALSE
  } else {
    has_convergence_issues <- length(messages[!grepl("singular", messages, ignore.case = TRUE)]) != 0
  }
  if (has_convergence_issues) {
    print(paste0(goi, " - Model didn't converge. Removed from the list."))
    model_tpm[[goi]] <- NULL
    
  } else {
    coeffs$tpm <- summary(model_tpm[[goi]])$coefficients %>%
      as.data.frame() %>%
      rownames_to_column("metric") %>%
      tibble() %>%
      mutate(gene = goi, .before = metric) %>%
      mutate(Item = metric, .before = metric) %>%
      mutate(singular = ifelse(isSingular(model_tpm[[goi]]), TRUE, FALSE)) %>%
      rbind(coeffs$tpm)
  }
}


# 
# test_tibble_presence <- list()
# model_presence <- list()
# coeffs$presence <- tibble()
# genes_of_interest <- c("Chitinase", "Glucosyltransferase", "PAPS reductase", "PnuC")
# for (goi in genes_of_interest) {
#   
#   test_tibble_presence[[goi]] <- gene_tpm %>%
#     filter(gene == goi) %>%
#     left_join(., metadata[c("Sample_ID", "Country", "Hive_ID", "Season", "Gut_part")], by = "Sample_ID") %>%
#     distinct() %>%
#     left_join(., nosema_relabund, by = "Sample_ID") %>%
#     mutate(Gut_part = factor(Gut_part, levels = c("rec", "ile", "mid"))) %>%
#     mutate(Vairimorpha_presence = ifelse(nosema_relabund > 0.01, 1 , 0),
#            log_tpm = log10(tpm)) %>%
#     filter(!is.infinite(log_tpm))
#   
#   model_presence[[goi]] <- glmer(Vairimorpha_presence ~ log_tpm + Gut_part + Season +
#                                    (1 | Hive_ID ), data = test_tibble_presence[[goi]],
#                                  family = binomial)
#   
#   has_convergence_issues <- FALSE
#   messages <- model_presence[[goi]]@optinfo$conv$lme4$messages
#   if (is.null(messages)) {
#     has_convergence_issues <- FALSE
#   } else {
#     has_convergence_issues <- length(messages[!grepl("singular", messages, ignore.case = TRUE)]) != 0
#   }
#   if (has_convergence_issues) {
#     print(paste0(goi, " - Simple model didn't converge. Removed from the list."))
#     model_presence[[goi]] <- NULL
#     
#   } else {
#     coeffs$presence <- summary(model_presence[[goi]])$coefficients %>%
#       as.data.frame() %>%
#       rownames_to_column("metric") %>%
#       tibble() %>%
#       mutate(gene = goi, .before = metric) %>%
#       mutate(Item = metric, .before = metric) %>%
#       mutate(singular = ifelse(isSingular(model_presence[[goi]]), TRUE, FALSE)) %>%
#       rbind(coeffs$presence)
#   }
# }
#
# No significant results, so dropped at this point


#####
# EXTRACT SLOPES
slopes <- list()
for (level in names(coeffs)) {
  
  temp_slope_tibble <- coeffs[[level]] %>%
    filter(metric %in% c("log_Vairimorpha_relabund")) %>%
    rename(raw_p_value = `Pr(>|t|)`) %>%
    mutate(raw_p_significant = case_when(raw_p_value <= 0.001 ~ "***",
                                         raw_p_value <= 0.01 ~ "**",
                                         raw_p_value <= 0.05 ~ "*",
                                         raw_p_value <= 0.075 ~ ".",
                                         .default = "n.s."
    )) %>%
    mutate(test_name = paste0(gene, "; ", Item), .before = gene)
  
  slopes[[level]] <- coeffs[[level]] %>%
    filter(metric == "(Intercept)") %>%
    mutate(test_name = paste0(gene, "; ", temp_slope_tibble$Item), .before = gene) %>%
    rename(intercept = Estimate,
           sd_intercept = `Std. Error`) %>%
    select(test_name, intercept, sd_intercept) %>%
    full_join(temp_slope_tibble, ., by = "test_name") %>%
    relocate(c(intercept, sd_intercept), .after = `Std. Error`)
  
}

##### 
# ADJUST P-VALUES

all_slopes <- bind_rows(slopes) %>%
  mutate(p_adjusted = p.adjust(raw_p_value, method = "BH")) %>%
  mutate(p_adjust_significant = case_when(p_adjusted <= 0.001 ~ "***",
                                          p_adjusted <= 0.01 ~ "**",
                                          p_adjusted <= 0.05 ~ "*",
                                          p_adjusted <= 0.075 ~ ".",
                                          .default = "n.s."))

#####
# MAKE PLOTS

# Significant results
lowest_highest <- nosema_relabund %>%
  filter(nosema_relabund > 0) %>%
  reframe(Item = "log_Vairimorpha_relabund",
          lowest = min(log10(nosema_relabund)),
          highest = max(log10(nosema_relabund)))


sig_tests <- all_slopes %>% 
  filter(p_adjusted < 0.05) %>%
  left_join(., lowest_highest, by = "Item") %>%
  mutate(effect = linear_effect_fun(s = Estimate, h = highest, l = lowest)) %>%
  group_by(gene) %>%
  mutate(y_stretching_factor = max(effect) / effect) %>%
  mutate(which_y_end_to_stretch = if_else(
    Estimate[y_stretching_factor == 1] > 0,
    "upper_end",
    "lower_end")) %>% 
  ungroup() %>%
  mutate(y_stretching_factor = case_when(
    which_y_end_to_stretch == "upper_end" & Estimate < 0 ~ 1/y_stretching_factor,
    which_y_end_to_stretch == "lower_end" & Estimate > 0 ~ 1/y_stretching_factor,
    .default = y_stretching_factor
  )) %>%
  select(-effect)

color_list <- list(dark = list(Chitinase = "#8B4513"),
                   bright = list(Chitinase = "#8B4513"))

tested_gene <- "Chitinase"
tested_item <- "log_Vairimorpha_relabund"

tpm_plot <- mixed_model_plot(filt_test_tibble = sig_tests,
                   transform_fun = linear_fun,
                   effect_fun = linear_effect_fun,
                   dark_col = color_list$dark[[tested_gene]],
                   bright_col = color_list$bright[[tested_gene]],
                   y_axis_label = "Log rel. gene abund.") +
  labs(x = "Log rel. *Vairimorpha* abund.") +
  theme(axis.title = element_markdown())

common_legend <- legend_factory(title = "Gene", 
                                items = names(color_list$dark),
                                colors = unlist(color_list$dark),
                                position = "bottom")

wrap_of_one <- tpm_plot / common_legend + plot_layout(heights = c(20,1))

# All results
all_tests_forest_plot <- all_slopes %>%
  mutate(axis_labels = fct_rev(fct_inorder(test_name))) %>%
  mutate(estimate = Estimate,
         error = `Std. Error`) %>%
  forest_plot(plot_title = "all tests")

#####
# DIAGNOSTICS

model_diagnostics <- list()
for (item in "log_Vairimorpha_relabund") {
  for (goi in names(model_tpm)) {
    slope_of_interest <- all_slopes %>%
      filter(Item == item,
             gene == goi) %>%
      reframe(test_name = test_name,
              is_significant = p_adjusted < 0.05)
    if (slope_of_interest$is_significant) {
      model_diagnostics[[item]][[goi]] <- diagnostics_linear_model(model_tpm[[goi]], name = slope_of_interest$test_name)
    }
  }
}

#####
# SAVE FILES

system("mkdir -p output/R/genes_pathogens_and_landuse/gene_tpm_vs_nosema_relabund/")

write_delim(coeffs$tpm, "output/R/genes_pathogens_and_landuse/gene_tpm_vs_nosema_relabund/gene_tpm_vs_nosema_relabund.all_coeffs.tsv",
            delim = "\t")
write_delim(all_slopes, "output/R/genes_pathogens_and_landuse/gene_tpm_vs_nosema_relabund/gene_tpm_vs_nosema_relabund.all_slopes.tsv",
            delim = "\t")
ggsave("output/R/genes_pathogens_and_landuse/gene_tpm_vs_nosema_relabund/gene_tpm_vs_nosema_relabund.all_tests.pdf",
       all_tests_forest_plot, width = 8, height = 4)

ggsave("output/R/genes_pathogens_and_landuse/gene_tpm_vs_nosema_relabund/gene_tpm_vs_nosema_relabund.wrap.pdf",
       wrap_of_one, width = 3.5, height = 3.5)

ggsave("output/R/genes_pathogens_and_landuse/gene_tpm_vs_nosema_relabund/single_panel.Citinase.Vairimorpha_relabund.pdf",
       tpm_plot, width = 3.5, height = 3.5)
ggsave("output/R/genes_pathogens_and_landuse/gene_tpm_vs_nosema_relabund/single_panel.common_legend.pdf",
       common_legend, width = 6, height = 1)

write_delim(nosema_relabund, "output/R/genes_pathogens_and_landuse/gene_tpm_vs_nosema_relabund/nosema_relabund.tsv",
            delim = "\t")

ggsave(paste0("output/R/genes_pathogens_and_landuse/gene_tpm_vs_nosema_relabund/model_diagnostics.log_Vairimorpha_relabund.Chitinase.pdf"),
       model_diagnostics$log_Vairimorpha_relabund$Chitinase, width = 6, height = 6)

saveRDS(tpm_plot, "output/R/genes_pathogens_and_landuse/gene_tpm_vs_nosema_relabund/single_panel.RDS.Citinase.Vairimorpha_relabund.rds")

