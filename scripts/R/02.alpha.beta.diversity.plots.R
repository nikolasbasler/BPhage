library(phyloseq)
library(furrr)
library(patchwork)
library(vegan)
library(ape)
library(ggpubr)
library(gridExtra)
library(VennDiagram)
library(tidyverse)
library(ggforce)
library(RColorBrewer)

start_time <- Sys.time()
source("scripts/R/helpers/tpm_and_reads_per_kb.R")
source("scripts/R/helpers/read_count_stats.R")
source("scripts/R/helpers/relative_abundance.R")
source("scripts/R/helpers/alpha_diversity.R")
source("scripts/R/helpers/beta_diversity.R")
source("scripts/R/helpers/decontam.R")
source("scripts/R/helpers/venn.R")

set.seed(1)
iterations <- 1000
# color_vector <- c(brewer.pal(n = 8, name = "Dark2"), brewer.pal(n = 12, name = "Set3"), brewer.pal(n = 8, name = "Set1")[c(1,2,8)])
# color_vector <- colorRampPalette(color_vector)(43)
color_vector <- colorRampPalette(brewer.pal(n = 8, name = "Set3"))(43) %>%
  sample()

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Read in data and generate abundance and TPM tables ####

metadata <- read.csv("data/metadata.csv") %>%
  mutate(Country = factor(Country, levels = c("PT", "FR", "UK", "BE", "NL", "CH", "DE", "RO"))) %>%
  mutate(Season = factor(Season, levels = c("spr", "sum", "aut"))) %>%
  mutate(Gut_part = factor(Gut_part, levels = c("mid", "ile", "rec"))) %>%
  mutate(Health = factor(Health, levels = c("very_bad", "bad", "good", "very_good")))
sample_order <- metadata %>% 
  arrange(Country, Season, Hive_ID, Gut_part) %>%
  select(Sample_ID) %>%
  unlist(use.names = FALSE)
metadata <- metadata %>%
  mutate(Sample_ID = factor(Sample_ID, levels=sample_order))
row.names(metadata) = metadata$Sample_ID
classification <- read.csv("output/R/euk_filtered/euk.diamond.and.blast.classification.csv")
contig_order <- classification %>%
  arrange(Phylum, Class, Order, Family, Subfamily, Genus, Species) %>% 
  select(contig) %>%
  unlist(use.names = FALSE)

classification <- classification %>%
  mutate(contig = factor(contig, levels=contig_order)) %>%
  filter(Species!="Apis mellifera filamentous virus")

euk_abundance <- read.csv("output/R/euk_filtered/euk.abundance.contig.csv") %>%
  select(!contains("Blank"))
blanks_abundance <- read.csv("output/R/euk_filtered/euk.abundance.contig.csv") %>%
  select(contig, contains("Blank"))


## Generate abundance tables for different taxonomic levels.####
taxlevels <- c("contig", "Species", "Genus", "Subfamily", "Family", "Order", "Class", "Phylum", "Host_group")
euk_ab <- list()
euk_ab = taxlevels %>%
  set_names() %>%
  map(~tax_sum(., ab_table=euk_abundance,
               classif = classification))

euk_lengths <- list()
euk_lengths <- taxlevels %>%
  set_names() %>%
  map(~tax_lengths(., classif = classification))

taxlevels <- c("contig", "Species", "Genus", "Family")
hostgroups <- classification$Host_group %>% unique()
euk_ab_hostgroup <- list()
for (hostg in hostgroups) {
  hostg_filt <- hostg_filter(hg=hostg, ab_table=euk_ab$contig,
          classif = classification)
  euk_ab_hostgroup[[hostg]] <- taxlevels %>%
    set_names() %>%
    map(~tax_sum(., ab_table=hostg_filt,
                 classif = classification))
}

## Generate TPM tables for different taxonomic levels. ####
euk_tpm <- list()
for (lvl in names(euk_ab)) {
  euk_tpm[[lvl]] <- euk_ab[[lvl]] %>% 
    inner_join(., euk_lengths[[lvl]], by=lvl) %>%
    calc_tpm(., lvl)
}
euk_tpm_hostgroup <- list()
for (hostg in hostgroups) {
  for (lvl in names(euk_ab_hostgroup[[hostg]])) {
    euk_tpm_hostgroup[[hostg]][[lvl]] <- euk_ab_hostgroup[[hostg]][[lvl]] %>%
      inner_join(., euk_lengths[[lvl]], by=lvl) %>% 
      calc_tpm(., lvl)
  }
}

## Make additional abundance tables with all read counts where the TPM table is below the threshold is set to 0. ####
tpm_thresh <- 0.01
euk_ab_filt <- list()
for (lvl in names(euk_ab)) {
  euk_ab_filt[[lvl]] <- euk_ab[[lvl]]
  if (!identical(rownames(euk_tpm[[lvl]]), rownames(euk_ab_filt[[lvl]])) |
      !identical(colnames(euk_tpm[[lvl]]), colnames(euk_ab_filt[[lvl]]))
      ) 
  {
    stop("Row or column names don't match up!")
  }
  euk_ab_filt[[lvl]][euk_tpm[[lvl]] < tpm_thresh] <- 0
}
euk_ab_hostgroup_filt <- list()
for (hostg in hostgroups) {
  for (lvl in names(euk_ab_hostgroup[[hostg]])) {
    euk_ab_hostgroup_filt[[hostg]][[lvl]] <- euk_ab_hostgroup[[hostg]][[lvl]]
    if (!identical(rownames(euk_tpm_hostgroup[[hostg]][[lvl]]), rownames(euk_ab_hostgroup_filt[[hostg]][[lvl]])) |
        !identical(colnames(euk_tpm_hostgroup[[hostg]][[lvl]]), colnames(euk_ab_hostgroup_filt[[hostg]][[lvl]]))
        ) 
    {
      stop("Row or column names don't match up!")
    }
    euk_ab_hostgroup_filt[[hostg]][[lvl]][euk_tpm_hostgroup[[hostg]][[lvl]] < tpm_thresh] <- 0
    
  }
}

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Cross-contamination plot from extraction plates using DWV ####

spec_tpm_with_blanks <- blanks_abundance %>%
  tax_sum("Species", ab_table=., classif = classification) %>%
  inner_join(., euk_lengths$Species, by="Species") %>%
  calc_tpm(., "Species") %>% 
  inner_join(., euk_tpm$Species, by="Species")

ccontam_tax <- "Deformed wing virus"
ccontam_p <- crosscontam_plot(spec_tpm_with_blanks, classification, ccontam_tax)

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Set min seq thresholds for rarefaction. ####

# count_stats_all <- report_stats(df = euk_ab$contig, 
#                             thresholds=c(5069, 5740)) # First number: quantile_90_ratio < 1000, second number: quantile_95_ratio < 2500
# min_seq_all <- 5740 # For quantile_90_ratio < 1000. 18 discarded samples (1 mid, 15 ile, 2 rec). If done on species level, the thresholds are exactly the same.
# discarded_on_all <- discards(count_stats_all$ratios, min_seq_all)$discarded
# lost_bees_on_all <- discards(count_stats_all$ratios, min_seq_all)$lost_bees # No bee pool lost completely. Only gut parts from differenct locations/time points.

count_stats_arthro <- report_stats(df = euk_ab_hostgroup$arthropod$contig,
                                 thresholds=c(8177, 8586)) # First number: quantile_90_ratio < 1000, second number: quantile_95_ratio < 2500
min_seq_arthro <- 8586 # For quantile_90_ratio < 1000. 34 discarded samples (2 mid 24 ile 8 rec). If done on species level, the thresholds are again exactly the same.
discarded_on_arthro <- discards(count_stats_arthro$ratios, min_seq_arthro)$discarded
lost_bees_on_arthro <- discards(count_stats_arthro$ratios, min_seq_arthro)$lost_bees # Again, no bee pool lost completely. Only gut parts from differenct locations/time points.

##
# 1% TPM filter: Same thrsholds as no TPM filter
# 2% TPM filter: Thresholds both 7752
#       min_seq_arthro_filt <- 8177 # For quantile_90_ratio < 1000. 35 discarded samples (3 mid 24 ile 8 rec). 
# 3% TPM filter: Thresholds both 7430
#       min_seq_arthro_filt <- 7752 # For quantile_90_ratio < 1000. 35 discarded samples (3 mid 25 ile 7 rec). 

# count_stats_arthro_filt <- report_stats(df = euk_ab_hostgroup_filt$arthropod$contig,
#                                    thresholds=c(7430, 7430)) # First number: quantile_90_ratio < 1000, second number: quantile_95_ratio < 2500
# min_seq_arthro_filt <- 7752 # For quantile_90_ratio < 1000. 35 discarded samples (3 mid 25 ile 7 rec). 
# discarded_on_arthro_filt <- discards(count_stats_arthro_filt$ratios, min_seq_arthro)$discarded
# lost_bees_on_arthro_filt <- discards(count_stats_arthro_filt$ratios, min_seq_arthro)$lost_bees # Here, one bee pool gets lost completely: BE_16562_spr



#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Alpha ####

met_v <- c("Country", "Season", "Gut_part", "Health")
taxlevels <- c("contig", "Species", "Genus", "Family")
alpha_arthro <- list()
for (tlvl in taxlevels) {
alpha_arthro[[tlvl]] <- alpha_stats(df = euk_ab_hostgroup_filt$arthropod[[tlvl]], 
                                  meta_vars = met_v, 
                                  min_seq = min_seq_arthro, # For 1% TPM filter
                                  df_lengths = euk_lengths[[tlvl]])
}

##
# richest <- alpha_arthro_species$table %>%
#   filter(Richness>=10)
# euk_tpm_hostgroup$arthropod$Species %>%
#   select(Species, all_of(richest$Sample_ID)) %>%
#   mutate(across(-Species, ~round(., digits=3))) %>% View()

##

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Beta ####

met_v <- c("Country", "Season", "Gut_part", "Health")
taxlevels <- c("contig", "Species", "Genus", "Family")
beta_dist_arthro <- list()
beta_plot_list_arthro <- list()
for (tlvl in taxlevels) {
  beta_dist_arthro[[tlvl]] <- rared_ordination(df = euk_ab_hostgroup_filt$arthropod[[tlvl]], 
                                      meta_vars = met_v, 
                                      min_seq = min_seq_arthro, # For 1% TPM filter
                                      df_lengths = euk_lengths[[tlvl]])
  
  beta_plot_list_arthro[[tlvl]] <- beta_plot(beta_dist_arthro[[tlvl]]$ord_list, met_v)
}

# beta_dist_arthro_contig <- rared_ordination(df = euk_ab_hostgroup_filt$arthropod$contig, 
#                                             meta_vars = met_v, 
#                                             min_seq = min_seq_arthro, # For 1% TPM filter
#                                             df_lengths = euk_lengths$contig)
# beta_dist_arthro_species <- rared_ordination(df = euk_ab_hostgroup_filt$arthropod$Species, 
#                                              meta_vars = met_v, 
#                                              min_seq = min_seq_arthro, # For 1% TPM filter
#                                              df_lengths = euk_lengths$Species)
# 
# 
# beta_plot_list_arthro_contig <- beta_plot(beta_dist_arthro_contig$ord_list, met_v)
# beta_plot_list_arthro_species <- beta_plot(beta_dist_arthro_species$ord_list, met_v)

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Contig overview heatmap ####

family_heatmaps_row <- list()
plotted_samples <- list()
plotted_contigs <- list()
for (fam in unique(classification$Family)) {
  if (fam == "") {
    next
  }
  # All countries in one row:
    fam_tpm <- euk_tpm$contig %>%
      inner_join(., classification, by="contig") %>% 
      filter(Family==fam) %>%
      select(colnames(euk_tpm$contig)) %>%
      select(which(colSums(. != 0) > 0))
    plotted_samples[[fam]] <- ncol(fam_tpm)
    plotted_contigs[[fam]] <- nrow(fam_tpm)
    family_heatmaps_row[[fam]] <- contig_heatmap(df = fam_tpm,
                                                   classif = classification)
}
split_families <- c("", "Parvoviridae", "Dicistroviridae")
for (fam in split_families) {
  species_names <- classification %>% 
    filter(Species!="Apis mellifera filamentous virus") %>%
    filter(Family==fam) %>%
    select(Species) %>%
    unlist(use.names=FALSE) %>%
    unique()
  fam_tpm <- euk_tpm$contig %>%
    inner_join(., classification, by="contig") %>% 
    filter(Family==fam) %>%
    select(colnames(euk_tpm$contig))
  if (fam=="") {
    fam <- "NoFamily"
  }
  for (speci in species_names) {
    spec_tpm <- fam_tpm %>%
      inner_join(., classification, by="contig") %>% 
      filter(Species==speci) %>%
      select(colnames(fam_tpm)) %>%
      select(which(colSums(. != 0) > 0))
    if (speci=="") {
      speci <- paste0(fam, "_NoSpecies")
    } else {
      speci <- paste0(fam, "_", speci)
    }
    plotted_samples[[speci]] <- ncol(spec_tpm)
    plotted_contigs[[speci]] <- nrow(spec_tpm)
    family_heatmaps_row[[speci]] <- contig_heatmap(df = spec_tpm,
                                                 classif = classification)
  }
}

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Host group histrograms: ####
hostgroup_hist <- list()
for (hostgroup in euk_tpm$Host_group$Host_group) {
  hostgroup_hist[[hostgroup]] <- euk_tpm$Host_group %>%
    # rename(group=name) %>%
    pivot_longer(-Host_group) %>%
    filter(Host_group==hostgroup) %>%
    ggplot(aes(x=value)) +
    geom_histogram(binwidth = 0.01) +
    ggtitle(paste0("TPMs of hostgroup ",hostgroup))
}

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Average TPM ####
met_v <- c("Country", "Season", "Gut_part", "Health", "Sample_ID")
tax_levels <- c("Species", "Genus", "Family")
average_tpm_plots <- list()
average_tpm_plots$Host_group <- average_tpm_bar_plot(tpm_table=euk_tpm$Host_group,
                        tl="Host_group", hg="all",
                        meta_vars=met_v)
for (hgr in names(euk_tpm_hostgroup)) {
  for (tlvl in tax_levels) {
    average_tpm_plots[[hgr]][[tlvl]] <- average_tpm_bar_plot(
      tpm_table = euk_tpm_hostgroup[[hgr]][[tlvl]],
      tl = tlvl,
      hg = hgr,
      meta_vars = met_v,
      threshold_for_other=0.03)
  }
}

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Prevalence ####
met_v <- c("Country", "Season", "Gut_part", "Health")
tax_levels <- c("Species", "Genus", "Family")
prevalence_plots <- list()
for (hgr in names(euk_ab_hostgroup_filt)) {
  for (tlvl in tax_levels) {
    trsh <- 3
    # if (tlvl=="Species") {
    #   trsh <- 10
    # } else {
    #   trsh <- 3
    # }
    prevalence_plots[[hgr]][[tlvl]] <- prevalence_bar_plot(
      abtable = euk_ab_hostgroup_filt[[hgr]][[tlvl]],
      tl = tlvl, 
      hg = hgr,
      meta_vars = met_v,
      threshold_for_other = trsh)
  }
}

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Country by country average TPM and prevalence: ####
tax_levels <- c("Species", "Genus", "Family")
average_tpm_plots_per_country <- list()
prevalence_plots_per_country <- list()
for (tlvl in tax_levels) {
  for (countr in levels(metadata$Country)) {
    tpm_filt <- euk_tpm_hostgroup$arthropod[[tlvl]] %>%
      select(all_of(tlvl), starts_with(countr))
    country_tpm_plots <- average_tpm_bar_plot(
      tpm_table = tpm_filt, 
      tl = tlvl, 
      hg = "arthropod", 
      meta_vars = c("Season", "Gut_part", "Health", "Sample_ID"),
      title_prefix = paste0(countr, " - "),
      threshold_for_other=0.03)
      average_tpm_plots_per_country[[tlvl]][[countr]] <- 
        country_tpm_plots$Sample_ID /
        (country_tpm_plots$Season + country_tpm_plots$Gut_part + country_tpm_plots$Health)
      
      ab_filt <- euk_ab_hostgroup_filt$arthropod[[tlvl]] %>%
          select(all_of(tlvl), starts_with(countr))
      number_of_samples <- ncol(ab_filt)-1
      prevalence_plots_per_country[[tlvl]][[countr]] <- prevalence_bar_plot(
        abtable = ab_filt,
        tl = tlvl, 
        hg = "arthropod",
        meta_vars = c("Season", "Gut_part", "Health"),
        title_prefix = paste0(countr, " (",number_of_samples," samples) - "),
        threshold_for_other=3)
  }
}

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Prevalenve Venn diagrams ####
# Output written to files already here.
venn_stats <- list()
tax_levels <- c("contig", "Species", "Genus", "Family")
met_v <- c("Season", "Gut_part", "Health")
for (tlvl in tax_levels) {
  for (m_var in met_v) {
    system(paste0("mkdir -p output/R/venns/",tlvl,"/",m_var))
    for (cntr in c("all", levels(metadata$Country))) {
      venn_stats[[tlvl]][[m_var]][[cntr]] <- prevalence_venn(
        abtable = euk_ab_hostgroup_filt$arthropod[[tlvl]],
        meta_var = m_var,
        country=cntr)
    }
  }
}

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Save files ####

## Cross-contamination DWV plot ####
ggsave(paste0("output/R/plate.cross.",ccontam_tax,".pdf"), ccontam_p, 
       width = 15, height = 10)

## Alpha diversity plots and tables ####
system("mkdir -p output/R/alpha")
for (tax in names(alpha_arthro)) {
  ggsave(paste0("output/R/alpha/alpha.arthro.",tax,".pdf"),
         alpha_arthro[[tax]]$plot, width = 12, height=10)
  write_csv(alpha_arthro[[tax]]$table,
            paste0("output/R/alpha/alpha.arthro.",tax,".csv"))
}

## Beta diversity plots and tables ####
for (tax in names(beta_plot_list_arthro)) {
  path <- paste0("output/R/beta/", tax, "_plots")
  system(paste0("mkdir -p ", path))
  for (p in names(beta_plot_list_arthro[[tax]])) {
    if (p=="all") {
      width <- 17
      height <- 15
    } else {
      width <- 17
      height <- 5
    }
    for (q in names(beta_plot_list_arthro[[tax]][[p]])) {
      ggsave(paste0(path,"/beta.arthro.", tax,".", p, ".", q, ".pdf"),
             beta_plot_list_arthro[[tax]][[p]][[q]], width = width, height=height)
    }
  }
  
  write_csv(rownames_to_column(beta_dist_arthro[[tax]]$avg_dist, "Sample_ID"),
            paste0("output/R/beta/beta_dist_arthro_", tax,".csv"))
}

## Heatmaps ####
system("mkdir -p output/R/heatmaps/")
for (tax in names(family_heatmaps_row)) {
  ggsave(paste0("output/R/heatmaps/", tax, ".pdf"),
         family_heatmaps_row[[tax]], 
         width = 10+plotted_samples[[tax]]/5, height =2+plotted_contigs[[tax]]/5, 
         limitsize = FALSE)
}

## TPM ####
system("mkdir -p output/R/relative_abundance_overall/")
for (p in names(hostgroup_hist)) {
  ggsave(paste0("output/R/relative_abundance_overall/hostgroup.hist.", p, ".pdf"),
         hostgroup_hist[[p]], width=8, height=4)
}
for (tl in names(average_tpm_plots$Host_group)) {
  ggsave(paste0("output/R/relative_abundance_overall/average.TPM.Host_groups.",tl,".pdf"),
         average_tpm_plots$Host_group[[tl]], width=15, height=8)
}

for (tl in names(average_tpm_plots$arthropod)) {
  for (me in names(average_tpm_plots$arthropod[[tl]])) {
    if (tl=="Species" && me=="Sample_ID") {
      wid=30
    } else {
      wid=15
    }
    ggsave(paste0("output/R/relative_abundance_overall/average.TPM.arthropod.",tl,".",me,".pdf"),
           average_tpm_plots$arthropod[[tl]][[me]], width=wid, height=8)
  }
}

## TPM per country ####
for (tl in names(average_tpm_plots_per_country)) {
  system(paste0("mkdir -p output/R/relative_abundance_per_country/",tl))
  if (tl=="Species") {
    wid=25
    hei=12
  } else {
    wid=15
    hei=12  
  }
  for (countr in names(average_tpm_plots_per_country[[tl]])) {
    ggsave(paste0("output/R/relative_abundance_per_country/", tl,"/average.TPM.country.",
                  countr,".",tl,".pdf"),
           average_tpm_plots_per_country[[tl]][[countr]], width=wid, height=hei)
  }
}

## Prevalence ####
system("mkdir -p output/R/prevalence_overall/")
for (tl in names(prevalence_plots$arthropod)) {
  if (tl=="Species") {
    wid=15
  } else {
    wid=10
  }
  ggsave(paste0("output/R/prevalence_overall/prevalence.arthropod.",tl,".pdf"),
         prevalence_plots$arthropod[[tl]], width=wid, height=20,
         limitsize = FALSE)
}

## Prevalence per country ####
for (tl in names(prevalence_plots_per_country)) {
  system(paste0("mkdir -p output/R/prevalence_per_country/",tl))
  if (tl=="Species") {
    wid=12
    hei=12
  } else {
    wid=8
    hei=12  
  }
  for (countr in names(prevalence_plots_per_country[[tl]])) {
    ggsave(paste0("output/R/prevalence_per_country/", tl,"/prevalence.arthropod.country.",
                  countr,".",tl,".pdf"),
           prevalence_plots_per_country[[tl]][[countr]], width=wid, height=hei)
  }
}
Sys.time() - start_time

