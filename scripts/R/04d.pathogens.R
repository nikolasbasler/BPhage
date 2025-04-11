library(readxl)
library(lme4)
library(lmerTest)
library(ggpubr)
library(ggrepel)
library(patchwork)
library(tidyverse)

source("scripts/R/helpers/mixed_helpers.R")

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
         across(starts_with("Cat"), ~ factor(.x, levels = c("L", "M", "H")))) %>%
  mutate(nosema_spores = ifelse(`N. spores` == "ND", "0", `N. spores`),
         nosema_spores = ifelse(`N. spores` == '< 25000', "25000", nosema_spores),
         nosema_spores = as.numeric(nosema_spores))

pathogens_Cts <- c("DWV A", "DWV B", "ABPV", "CBPV", "BQCV", "SBV", "EFB",
               "AFB", "AFB (cfu)", "N. apis", "N. ceranae")

pathogens_cats <- c("Cat DWV A", "Cat DWV B", "Cat ABPV", "Cat CBPV", "Cat BQCV",
               "Cat SBV", "Cat EFB", "Cat AFB")

meta_with_pathogens <- pathogen_data %>% 
  filter(Years == 2020) %>%
  select(BGOOD_sample_code, all_of(c(pathogens_Cts, pathogens_cats, "nosema_spores"))) %>%
  left_join(metadata[c("Sample_ID", "BGOOD_sample_code", "Bee_pool", "Country", "Hive_ID", "Season", "Gut_part")], ., by = "BGOOD_sample_code") %>%
  distinct() %>%
  mutate(across(all_of(pathogens_Cts), ~ str_replace_all(., "negative", "40"))) %>%
  mutate(across(all_of(pathogens_Cts), ~ suppressWarnings(as.numeric(.))))

mean_gene_tpm_and_pathogens <- gene_tpm %>% 
  left_join(., meta_with_pathogens, by = "Sample_ID") %>%
  # group_by(gene, Bee_pool) %>%
  # mutate(mean_tpm = mean(tpm), .before = BGOOD_sample_code) %>% 
  # ungroup() %>%
  select(gene, Sample_ID, Bee_pool, Country, Hive_ID, Season, Gut_part, tpm, all_of(c(pathogens_Cts, pathogens_cats, "nosema_spores"))) %>%
  distinct()

informative_numbers <- tibble()
for (pathogen in pathogens_Cts) {
  NAs <- mean_gene_tpm_and_pathogens %>% 
    filter(is.na(.data[[pathogen]])) %>%
    nrow()
  zeros <- mean_gene_tpm_and_pathogens %>%
    filter(.data[[pathogen]] == 40) %>%
    nrow()
  informative_numbers <- tibble(pat = pathogen, 
         NA_prop = NAs / nrow(mean_gene_tpm_and_pathogens), 
         zero_prop = zeros / nrow(mean_gene_tpm_and_pathogens),
         informative_numbers_prop = (nrow(mean_gene_tpm_and_pathogens) - NAs - zeros) / nrow(mean_gene_tpm_and_pathogens)) %>%
    rbind(informative_numbers, .)
}

NAs <- mean_gene_tpm_and_pathogens %>% 
  filter(is.na(nosema_spores)) %>%
  nrow()
zeros <- mean_gene_tpm_and_pathogens %>%
  filter(nosema_spores == 0) %>%
  nrow()
informative_numbers <- tibble(pat = "nosema_spores", 
       NA_prop = NAs / nrow(mean_gene_tpm_and_pathogens), 
       zero_prop = zeros / nrow(mean_gene_tpm_and_pathogens),
       informative_numbers_prop = (nrow(mean_gene_tpm_and_pathogens) - NAs - zeros) / nrow(mean_gene_tpm_and_pathogens)) %>%
  rbind(informative_numbers, .)

mean_gene_tpm_and_pathogens %>% 
  select(Bee_pool, nosema_spores) %>%
  distinct() %>%
  filter(!is.na(nosema_spores)) %>%
  ggplot(aes(x = nosema_spores)) +
  geom_histogram(bins = 50)

informative_numbers %>%
  arrange(desc(informative_numbers_prop)) %>%
  filter(!str_detect(pat, "Cat"))

gene_pathogen_tpm_tibble <- list()
model_tpm <- list()
coeffs_tpm <- tibble()

#####
# TPM TEST
pathogens_of_interest <- c("BQCV", "SBV", "DWV B")
# for (goi in unique(gene_tpm$gene)) {
for (goi in "phosphoadenosine phosphosulfate reductase") {
  for (poi in pathogens_of_interest) {
    # goi = "phosphoadenosine phosphosulfate reductase"
    # poi = "BQCV"
    
    gene_pathogen_tpm_tibble[[goi]][[poi]] <- mean_gene_tpm_and_pathogens %>%
      filter(gene == goi) %>%
      rename(pathogen_ct = all_of(poi)) %>% 
      filter(!is.na(pathogen_ct)) %>%
      select(Sample_ID, Bee_pool, Country, Hive_ID, Season, Gut_part, tpm, pathogen_ct) %>%
      # mutate(logged_ct = log2(pathogen_ct)) %>%
      mutate(log_tpm = log10(tpm)) %>%
      filter(!is.infinite(log_tpm))
    
    model_tpm[[goi]][[poi]] <- lmer(log_tpm ~ pathogen_ct + Season + Gut_part + 
           (1 | Bee_pool + Hive_ID ),
         data = gene_pathogen_tpm_tibble[[goi]][[poi]])
    
    summary(model_tpm[[goi]][[poi]])
    
    coeffs_tpm <- summary(model_tpm[[goi]][[poi]])$coefficients %>%
      as.data.frame() %>%
      rownames_to_column("metric") %>%
      tibble() %>%
      mutate(gene = goi, 
             pathogen = poi,
             .before = metric) %>%
      mutate(singular = ifelse(isSingular(model_tpm[[goi]][[poi]]), TRUE, FALSE)) %>%
      rbind(coeffs_tpm)
  }
}

summary(model_tpm$`phosphoadenosine phosphosulfate reductase`$BQCV)
summary(model_tpm$`phosphoadenosine phosphosulfate reductase`$SBV)
summary(model_tpm$`phosphoadenosine phosphosulfate reductase`$`DWV B`)

plot(model_tpm$`phosphoadenosine phosphosulfate reductase`$BQCV, which = 1)
qqnorm(resid(model_tpm$`phosphoadenosine phosphosulfate reductase`$BQCV))

lowest_highest_tpm <- mean_gene_tpm_and_pathogens %>%
  select(Sample_ID, all_of(pathogens_Cts)) %>%
  distinct() %>%
  pivot_longer(-Sample_ID, names_to = "pathogen", values_to = "ct") %>%
  filter(!is.na(ct)) %>%
  group_by(pathogen) %>%
  summarise(lowest = min(ct),
            highest = max(ct))

slopes_tpm <- coeffs_tpm %>%
  filter(metric == "pathogen_ct") %>%
  mutate(p_adjusted = p.adjust(`Pr(>|t|)`, method = "BH")) %>%
  mutate(p_adj_sig = case_when(p_adjusted <= 0.001 ~ "***",
                               p_adjusted <= 0.01 ~ "**",
                               p_adjusted <= 0.05 ~ "*",
                               p_adjusted <= 0.075 ~ ".",
                               .default = "n.s."
  )) %>%
  left_join(., lowest_highest_tpm, by = "pathogen") %>%
  mutate(change_in_range = Estimate * (highest - lowest),
         fold_change_in_range = 10^change_in_range,
         backtrans_estimate = 10^Estimate-1,
         backtrans_error = 10^Estimate - 10^(Estimate - `Std. Error`)) %>%
  mutate(axis_labels = pathogen) %>%
  mutate(axis_labels = fct_rev(fct_inorder(axis_labels)))

  

lme_plot <- forest_plot(slopes_tpm, plot_title = "Sulf phage relative abundanve ~ pathogen Ct")

fold_change_in_range_plot <- slopes_tpm %>%
  ggplot(aes(x = axis_labels, y = fold_change_in_range)) +
  geom_col() +
  coord_flip() +
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        axis.text.y=element_blank())

pathogen_wrap <- lme_plot + fold_change_in_range_plot
pathogen_wrap

# 
# #####
# # LOGIT TEST (nothing is significant but maybe try to correlate this rather to land use, not to tpm of sulf phages)
# pathogens_of_interest <- c("BQCV", "SBV", "DWV B",
#                            "ABPV", "N. ceranae", "DWV A", "CBPV")
# 
# gene_pathogen_presence_tibble <- list()
# model_logit <- list()
# coeffs_logit <- list()
# 
# for (goi in "phosphoadenosine phosphosulfate reductase") {
#   for (poi in pathogens_of_interest) {
#     # goi = "phosphoadenosine phosphosulfate reductase"
#     # poi = "ABPV"
#     
#     gene_pathogen_presence_tibble[[goi]][[poi]] <- mean_gene_tpm_and_pathogens %>%
#       filter(gene == goi) %>%
#       rename(pathogen_presence = all_of(poi)) %>% 
#       filter(!is.na(pathogen_presence)) %>%
#       mutate(pathogen_presence = ifelse(pathogen_presence < 40, 1, 0)) %>%
#       select(Sample_ID, Bee_pool, Country, Hive_ID, Season, Gut_part, tpm, pathogen_presence) %>%
#       mutate(log_tpm = log10(tpm)) %>%
#       filter(!is.infinite(log_tpm))
#     
#     # model_logit[[goi]][[poi]] <- glmer(pathogen_presence ~ log_tpm + Season + Gut_part +
#     #                               (1 | Bee_pool + Hive_ID ),
#     #                             data = gene_pathogen_presence_tibble[[goi]][[poi]],
#     #                             family = binomial)
#     
#     model_logit[[goi]][[poi]] <- lmer(log_tpm ~ pathogen_presence + Season + Gut_part +
#                                          (1 | Bee_pool + Hive_ID ),
#                                        data = gene_pathogen_presence_tibble[[goi]][[poi]])
# 
#     # summary(model_logit[[goi]][[poi]])
#     
#     coeffs_logit <- summary(model_logit[[goi]][[poi]])$coefficients %>%
#       as.data.frame() %>%
#       rownames_to_column("metric") %>%
#       tibble() %>%
#       mutate(gene = goi, 
#              pathogen = poi,
#              .before = metric) %>%
#       mutate(singular = ifelse(isSingular(model_logit[[goi]][[poi]]), TRUE, FALSE)) %>%
#       rbind(coeffs_logit)
#   }
# }
# 
# 
# slopes_logit <- coeffs_logit %>%
#   filter(metric == "pathogen_presence") %>%
#   mutate(p_adjusted = p.adjust(`Pr(>|t|)`, method = "BH")) %>%
#   mutate(p_adj_sig = case_when(p_adjusted <= 0.001 ~ "***",
#                                p_adjusted <= 0.01 ~ "**",
#                                p_adjusted <= 0.05 ~ "*",
#                                p_adjusted <= 0.075 ~ ".",
#                                .default = "n.s."
#   ))
# 
# slopes_logit <- coeffs_logit %>%
#   filter(metric == "log_tpm") %>%
#   mutate(p_adjusted = p.adjust(`Pr(>|z|)`, method = "BH")) %>%
#   mutate(p_adj_sig = case_when(p_adjusted <= 0.001 ~ "***",
#                                p_adjusted <= 0.01 ~ "**",
#                                p_adjusted <= 0.05 ~ "*",
#                                p_adjusted <= 0.075 ~ ".",
#                                .default = "n.s."
#   ))
# 
#   # left_join(., lowest_highest_logit, by = "pathogen") %>%
#   # mutate(change_in_range = Estimate * (highest - lowest),
#   #        fold_change_in_range = 10^change_in_range,
#   #        backtrans_estimate = 10^Estimate-1,
#   #        backtrans_error = 10^Estimate - 10^(Estimate - `Std. Error`)) %>%
#   # mutate(axis_labels = pathogen) %>%
#   # mutate(axis_labels = fct_rev(fct_inorder(axis_labels)))
# 


#####
# SAVE FILES
system("mkdir -p output/R/gene_content/pathogens")
ggsave("output/R/gene_content/pathogens/pathogen_wrap.pdf",
       pathogen_wrap, width = 7, height = 3)





#
# 
# pathogens_of_interest <- c("BQCV", "SBV", "DWV B")
# # pathogens_of_interest <- c("BQCV", "SBV", "DWV B", "ABPV", "N. ceranae", "DWV A", "CBPV", "N. spores")
# gene_pathogen_cor_tibble <- list()
# gene_pathogen_cor_plot <- list()
# gene_pathogen_cor_test <- list()
# gene_pathogen_cor_p_values <- tibble()
# # for (goi in unique(gene_tpm$gene)) {
# for (goi in "phosphoadenosine phosphosulfate reductase") {
#     for (poi in pathogens_of_interest) {
#       # goi = "phosphoadenosine phosphosulfate reductase"
#       # poi = "BQCV"
#     gene_pathogen_cor_tibble[[goi]][[poi]] <- mean_gene_tpm_and_pathogens %>%
#       filter(gene == goi) %>%
#       filter(!is.na(.data[[poi]])) %>%
#       select(Sample_ID, Country, Hive_ID, Season, tpm, all_of(poi)) %>%
#       rename(pathogen_ct = all_of(poi)) %>%
#       # mutate(logged_ct = log2(.data[[poi]])) %>%
#       mutate(log_tpm = log10(tpm)) %>%
#       filter(!is.infinite(log_tpm))
#     
#     if (nrow(gene_pathogen_cor_tibble[[goi]][[poi]]) > 2) {
#       gene_pathogen_cor_plot[[goi]][[poi]] <- gene_pathogen_cor_tibble[[goi]][[poi]] %>%
#         # ggplot(aes(x = .data[[poi]], y = log_mean_tpm, color = Country)) +
#         # ggplot(aes(x = .data[[poi]], y = log_mean_tpm)) +
#         ggplot(aes(x = pathogen_ct, y = log_tpm)) +
#         # ggplot(aes(x = .data[[poi]], y = mean_tpm)) +
#         geom_point() +
#         geom_smooth(method = "glm", formula = y ~ x) +
#         # stat_cor(method = "spearman") +
#         stat_cor(method = "pearson") +
#         ggtitle(paste0(goi, " - ", poi))
#       gene_pathogen_cor_test[[goi]][[poi]] <- cor.test(gene_pathogen_cor_tibble[[goi]][[poi]]$pathogen_ct,  
#                gene_pathogen_cor_tibble[[goi]][[poi]]$log_tpm,
#                # method = "spearman")
#                method = "pearson")
# gene_pathogen_cor_p_values <- tibble(
#         test = paste0(goi, " - ", poi),
#         R = gene_pathogen_cor_test[[goi]][[poi]]$estimate,
#         p_value = gene_pathogen_cor_test[[goi]][[poi]]$p.value) %>%
#         rbind(gene_pathogen_cor_p_values, .)
#     }
#   }
# }
# 
# wrap_2 <- wrap_plots(gene_pathogen_cor_plot$`phosphoadenosine phosphosulfate reductase`, ncol =2)
# ggsave("output/R/gene_content/pathogens/pathogen_wrap_2.pdf",
#        wrap_2, width = 10, height = 8)
# 
# gene_pathogen_cor_p_values %>%
#   mutate(adjusted_p = p.adjust(p_value, method = "BH")) %>%
#   arrange(adjusted_p)
