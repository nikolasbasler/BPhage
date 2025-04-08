library(lme4)
library(lmerTest)
library(tidyverse)

metadata <- readRDS("output/R/R_variables/metadata.RDS") %>% 
  tibble() %>%
  mutate(Hive_ID = as.character(Hive_ID))

gene_tpm <- read.delim("output/R/gene_content/landuse/gene_tpm.tsv") %>% tibble()

nosema_mapped_counts <- read.delim("output/nosema_mapped_counts_all.tsv") %>% tibble()

test_tibble_tpm <- nosema_mapped_counts %>%
  mutate(nosema_relabund = mapped_to_Nosema / Hostout_R1_plus_R2) %>%
  mutate(nosema_log_relabund = log10(nosema_relabund)) %>%
  left_join(., gene_tpm, by = "Sample_ID") %>%
  filter(gene == "phosphoadenosine phosphosulfate reductase") %>%
  mutate(gene_log_tpm = log10(tpm)) %>%
  filter(!is.infinite(nosema_log_relabund),
         !is.infinite(gene_log_tpm)) %>%
  left_join(., metadata[c("Sample_ID", "Country", "Hive_ID", "Season", "Gut_part")], by = "Sample_ID") %>%
  select(Sample_ID, Country, Hive_ID, Season, Gut_part, nosema_log_relabund, gene_log_tpm) 

# test_result_tpm <- lmer(gene_log_tpm ~ nosema_log_relabund + Gut_part + Season +
#        (1 | Hive_ID), data = test_tibble_tpm)
test_result_tpm <- lmer(gene_log_tpm ~ nosema_log_relabund + Gut_part + Season +
                          (1 | Hive_ID ), data = test_tibble_tpm)
summary(test_result_tpm)

plot(test_result_tpm, which = 1)
qqnorm(resid(test_result_tpm))



# 
# test_tibble_presence <- nosema_mapped_counts %>%
#   mutate(nosema_presence = ifelse(mapped_to_Nosema > 100, 1, 0)) %>%
#   left_join(., gene_tpm, by = "Sample_ID") %>%
#   filter(gene == "phosphoadenosine phosphosulfate reductase") %>%
#   mutate(gene_log_tpm = log10(tpm)) %>%
#   filter(!is.infinite(gene_log_tpm)) %>%
#   left_join(., metadata[c("Sample_ID", "Country", "Hive_ID", "Season", "Gut_part")], by = "Sample_ID") %>%
#   select(Sample_ID, Country, Hive_ID, Season, Gut_part, nosema_presence, gene_log_tpm)
# 
# test_result_logit <- glmer(nosema_presence ~ gene_log_tpm + Gut_part + Season +
#                           (1 | Hive_ID), data = test_tibble_presence,
#                           family = binomial)
# summary(test_result_logit)

