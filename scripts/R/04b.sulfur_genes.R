library(tidyverse)
library(ggpubr)
library(multcompView)

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

# Nosema
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
