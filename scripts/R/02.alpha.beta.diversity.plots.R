start_time <- Sys.time()
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
row.names(metadata) <- metadata$Sample_ID
classification <- read.csv("output/R/phage.filt.gnmd.classification.csv") %>%
  mutate(contig_length = contig_length/1000) %>%
  rename(length_kb = contig_length)
contig_order <- classification %>%
  arrange(Realm, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>% 
  select(contig) %>%
  unlist(use.names = FALSE)

classification <- classification %>%
  mutate(contig = factor(contig, levels=contig_order))

phage_abundance <- read.csv("output/R/phage.filt.abundance.contig.csv") %>%
  select(!contains("Blank"))

# Make a pie chart of the different taxons
taxon_pie <- table(classification$lowest_taxon) %>%
  as.data.frame() %>%
  rename(Lowest_taxon = Var1, Genomes = Freq) %>%
  arrange(desc(Genomes)) %>%
  mutate(Lowest_taxon = factor(Lowest_taxon, levels = Lowest_taxon)) %>%
  ggplot(aes(x = "", y=Genomes, fill = Lowest_taxon)) +
  geom_bar(stat = "identity", color="black", linewidth=0.2) +
  coord_polar("y") +
  theme_classic() +
  labs(x = NULL, y = NULL) +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  scale_fill_brewer(palette="Set3")

## Generate abundance tables for different taxonomic levels.####
taxlevels <- c("contig", "Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom", "Realm", "Host_group")
phage_ab <- list()
phage_ab <- taxlevels %>%
  set_names() %>%
  map(~tax_sum(., ab_table=phage_abundance,
               classif = classification))
phage_lengths <- list()
phage_lengths <- taxlevels %>%
  set_names() %>%
  map(~tax_lengths(., classif = classification))

taxlevels <- c("contig", "Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom", "Realm")
hostgroups <- classification$Host_group %>% unique()
phage_ab_hostgroup <- list()
for (hostg in hostgroups) {
  hostg_filt <- hostg_filter(hg = hostg, ab_table = phage_ab$contig,
          classif = classification)
  phage_ab_hostgroup[[hostg]] <- taxlevels %>%
    set_names() %>%
    map(~tax_sum(., ab_table=hostg_filt,
                 classif = classification))
}

## Generate contig abundance table with merged metadata variables 
# (e.g. all gut parts or all samples from the same country merged). Used 
# for prealenve plots
meta_merges <- list(Bee_pools = "Gut_part", 
                    Countries = c("Season", "Gut_part", "Hive_ID"),
                    Hives = c("Season", "Gut_part"))
phage_ab_meta_merges <- list()
for (merge in names(meta_merges)) {
  phage_ab_meta_merges[[merge]] <- meta_sum(ab_table = phage_ab$contig, 
                                            meta_vars = meta_merges[[merge]])
}

## Generate TPM tables for different taxonomic levels. ####
phage_tpm <- list()
for (lvl in names(phage_ab)) {
  phage_tpm[[lvl]] <- calc_tpm(abtable = phage_ab[[lvl]], level = lvl, lengths_df = phage_lengths[[lvl]])
}

phage_tpm_hostgroup <- list()
for (hostg in hostgroups) {
  for (lvl in names(phage_ab_hostgroup[[hostg]])) {
    phage_tpm_hostgroup[[hostg]][[lvl]] <- calc_tpm(abtable = phage_ab_hostgroup[[hostg]][[lvl]], level = lvl, lengths_df = phage_lengths[[lvl]])
  }
}

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Set min seq thresholds for rarefaction. ####

# 50% horizontal coverage cutoff
# count_stats <- report_stats(df = phage_abundance,
#                             thresholds=c(2689, 1002)) # First number: quantile_90_ratio < 1000, second number: quantile_95_ratio < 2500
# min_seq_count <- 2713 # For quantile_90_ratio < 1000. 11 discarded samples (10 mid, 1 ile).

# 70% horizontal coverage cutoff
count_stats <- report_stats(df = phage_abundance,
                            thresholds=c(2437, 522)) # First number: quantile_90_ratio < 1000, second number: quantile_95_ratio < 2500
count_stats
min_seq_count <- 2442 # For quantile_90_ratio < 1000. 12 discarded samples (9 mid, 3 ile).
discarded <- discards(count_stats$ratios, min_seq_count)$discarded
lost_bees <- discards(count_stats$ratios, min_seq_count)$lost_bees # No bee pool lost completely. Only gut parts from different locations/time points.

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Alpha ####
met_v <- c("Country", "Season", "Gut_part", "Health")
# taxlevels <- c("contig", "Species", "Genus", "Family")
taxlevels <- c("contig")
alpha <- list()
for (tlvl in taxlevels) {
  alpha[[tlvl]] <- alpha_stats(df = phage_ab[[tlvl]], 
                               meta_vars = met_v, 
                               min_seq = min_seq_count,
                               df_lengths = phage_lengths[[tlvl]])
}

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Beta ####
met_v <- c("Country", "Season", "Gut_part", "Health")
# taxlevels <- c("contig", "Species", "Genus", "Family")
taxlevels <- c("contig")
beta_dist <- list()
beta_plot_list <- list()
dist_hist_list <- list()
for (tlvl in taxlevels) {
  beta_dist[[tlvl]] <- rared_ordination(df = phage_ab[[tlvl]], 
                                        meta_vars = met_v, 
                                        min_seq = min_seq_count,
                                        df_lengths = phage_lengths[[tlvl]])
  
  beta_plot_list[[tlvl]] <- beta_plot(ordination_list = beta_dist[[tlvl]]$ord_list, 
                                      meta_vars = met_v)
}

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
    fam_tpm <- phage_tpm$contig %>%
      inner_join(., classification, by="contig") %>%
      filter(Family==fam) %>%
      select(colnames(phage_tpm$contig)) %>%
      select(which(colSums(. != 0) > 0))
    plotted_samples[[fam]] <- ncol(fam_tpm)
    plotted_contigs[[fam]] <- nrow(fam_tpm)
    family_heatmaps_row[[fam]] <- contig_heatmap(df = fam_tpm,
                                                   classif = classification)
}
# split_families <- c("", "Parvoviridae", "Dicistroviridae")
# for (fam in split_families) {
#   species_names <- classification %>%
#     filter(Species!="Apis mellifera filamentous virus") %>%
#     filter(Family==fam) %>%
#     select(Species) %>%
#     unlist(use.names=FALSE) %>%
#     unique()
#   fam_tpm <- euk_tpm$contig %>%
#     inner_join(., classification, by="contig") %>%
#     filter(Family==fam) %>%
#     select(colnames(euk_tpm$contig))
#   if (fam=="") {
#     fam <- "NoFamily"
#   }
#   for (speci in species_names) {
#     spec_tpm <- fam_tpm %>%
#       inner_join(., classification, by="contig") %>%
#       filter(Species==speci) %>%
#       select(colnames(fam_tpm)) %>%
#       select(which(colSums(. != 0) > 0))
#     if (speci=="") {
#       speci <- paste0(fam, "_NoSpecies")
#     } else {
#       speci <- paste0(fam, "_", speci)
#     }
#     plotted_samples[[speci]] <- ncol(spec_tpm)
#     plotted_contigs[[speci]] <- nrow(spec_tpm)
#     family_heatmaps_row[[speci]] <- contig_heatmap(df = spec_tpm,
#                                                  classif = classification)
#   }
# }

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Host group histrograms: ####
hostgroup_hist <- list()
for (hostgroup in phage_tpm$Host_group$Host_group) {
  hostgroup_hist[[hostgroup]] <- phage_tpm$Host_group %>%
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
tax_levels <- c("contig", "Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom", "Realm")
average_tpm_plots <- list()
average_tpm_plots$Host_group <- average_tpm_bar_plot(tpm_table = phage_tpm$Host_group,
                        tl="Host_group", hg="all",
                        meta_vars=met_v)
for (hgr in names(phage_tpm_hostgroup)) {
  for (tlvl in tax_levels) {
    average_tpm_plots[[hgr]][[tlvl]] <- average_tpm_bar_plot(
      tpm_table = phage_tpm_hostgroup[[hgr]][[tlvl]],
      tl = tlvl,
      hg = hgr,
      meta_vars = met_v,
      threshold_for_other=0.01)
  }
}

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Prevalence ####
prevalence_histo <- list()
for (merge in names(phage_ab_meta_merges)) {
  prevalence_histo[[merge]] <- prevalence_histogram(abtable = phage_ab_meta_merges[[merge]],
                                                          plot_title = merge)
}

# 
# met_v <- c("Country", "Season", "Gut_part", "Health")
# tax_levels <- c("contig", "Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom", "Realm")
# prevalence_plots <- list()
# for (hgr in names(phage_tpm_hostgroup)) {
#   for (tlvl in tax_levels) {
#     # trsh <- 3
#     if (tlvl=="contig") {
#       trsh <- 20
#     } else {
#       trsh <- 1
#     }
#     prevalence_plots[[hgr]][[tlvl]] <- prevalence_bar_plot(
#       abtable = phage_tpm_hostgroup[[hgr]][[tlvl]],
#       tl = tlvl, 
#       hg = hgr,
#       meta_vars = met_v,
#       threshold_for_other = trsh)
#   }
# }



#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Country by country average TPM and prevalence: ####
tax_levels <- c("Species", "Genus", "Family")
average_tpm_plots_per_country <- list()
prevalence_plots_per_country <- list()
for (tlvl in tax_levels) {
  for (countr in levels(metadata$Country)) {
    tpm_filt <- phage_tpm_hostgroup$all[[tlvl]] %>%
      select(all_of(tlvl), starts_with(countr))
    country_tpm_plots <- average_tpm_bar_plot(
      tpm_table = tpm_filt, 
      tl = tlvl, 
      hg = "all", 
      meta_vars = c("Season", "Gut_part", "Health", "Sample_ID"),
      title_prefix = paste0(countr, " - "),
      threshold_for_other=0.03)
      average_tpm_plots_per_country[[tlvl]][[countr]] <- 
        country_tpm_plots$Sample_ID /
        (country_tpm_plots$Season + country_tpm_plots$Gut_part + country_tpm_plots$Health)
      
      ab_filt <- phage_ab_hostgroup$all[[tlvl]] %>%
          select(all_of(tlvl), starts_with(countr))
      number_of_samples <- ncol(ab_filt)-1
      prevalence_plots_per_country[[tlvl]][[countr]] <- prevalence_bar_plot(
        abtable = ab_filt,
        tl = tlvl, 
        hg = "all",
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
        abtable = phage_ab_hostgroup$all[[tlvl]],
        meta_var = m_var,
        country=cntr)
    }
  }
}

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# Save files ####

## Taxon pie
ggsave("output/R/taxon_pie.pdf", taxon_pie, width = 7, height = 5)

## Alpha diversity plots and tables ####
system("mkdir -p output/R/alpha")
for (tax in names(alpha)) {
  ggsave(paste0("output/R/alpha/alpha.",tax,".pdf"),
         alpha[[tax]]$plot, width = 12, height=10)
  write_csv(alpha[[tax]]$table,
            paste0("output/R/alpha/alpha.",tax,".csv"))
}

## Beta diversity plots and tables ####
for (tax in names(beta_plot_list)) {
  pcoa_path <- paste0("output/R/beta/", tax, "_pcoa")
  hist_path <- paste0("output/R/beta/", tax, "_hist")
  system(paste0("mkdir -p ", pcoa_path))
  system(paste0("mkdir -p ", hist_path))
  for (p in names(beta_plot_list[[tax]])) {
    if (p=="all") {
      width <- 17
      height <- 15
    } else {
      width <- 17
      height <- 5
    }
    for (q in names(beta_plot_list[[tax]][[p]])) {
      ggsave(paste0(pcoa_path,"/beta.", tax,".", p, ".", q, ".pcoa.pdf"),
             beta_plot_list[[tax]][[p]][[q]], width = width, height = height)
      
      if (q == "control") { next }
      ggsave(paste0(hist_path,"/beta.", tax,".", p, ".", q, ".hist.pdf"),
             beta_dist[[tax]]$dist_hist_list[[p]][[q]], width = 8, height = 6)
    }
  }
  
  write_csv(rownames_to_column(beta_dist[[tax]]$avg_dist, "Sample_ID"),
            paste0("output/R/beta/beta_dist_", tax,".csv"))
}

# ## Heatmaps ####
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

for (tl in names(average_tpm_plots$all)) {
  for (me in names(average_tpm_plots$all[[tl]])) {
    if (tl=="Species" && me=="Sample_ID") {
      wid=30
    } else {
      wid=15
    }
    ggsave(paste0("output/R/relative_abundance_overall/average.TPM.all.",tl,".",me,".pdf"),
           average_tpm_plots$all[[tl]][[me]], width=wid, height=8)
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
system("mkdir -p output/R/prevalence/")
for (tl in names(prevalence_histo)) {
  if (tl=="Bee_pools") {
    wid=20
    plot <- prevalence_histo[[tl]]$plot # + scale_x_continuous(breaks = seq(3,150,3))
  } else if (tl=="Hives") {
    wid=12
    plot <- prevalence_histo[[tl]]$plot 
  } else {
    wid=8
    plot <- prevalence_histo[[tl]]$plot
  }
  ggsave(paste0("output/R/prevalence/prevalence.",tl,".pdf"),
         plot, width=wid, height=5)
  write_csv(prevalence_histo[[tl]]$table, paste0("output/R/prevalence/prevalence.",tl,".csv"))
}

end_time <- Sys.time()
end_time - start_time

