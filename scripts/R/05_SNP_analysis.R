library(ape)
library(patchwork)
library(ggtree)
library(TreeDist)
library(tidyverse)

set.seed(1)

metadata <- readRDS("output/R/R_variables/metadata.RDS")
# datasets <- c("Acer_CO1", "Acer_genome", "Acer_mt", "host", "phages", "MAG_mapped", "core", "all_terLs", "core_terLs")
datasets <- c("Acer_genome", "host", "all_terLs", "phages", "16S", "MAG_mapped")

ska_distance <- list()
for (dataset in datasets) {
  ska_distance[[dataset]] <- read.delim(paste0("output/SNP_analysis/", dataset, ".distances.tsv")) %>%
    separate_wider_delim(Sample.1, "_", names = c("country", "hive", "season"), too_many = "drop") %>%
    mutate(Sample1 = paste(country, hive, season, sep="_")) %>%
    select(-country, -hive, -season) %>%
    separate_wider_delim(Sample.2, "_", names = c("country", "hive", "season"), too_many = "drop") %>%
    mutate(Sample2 = paste(country, hive, season, sep="_")) %>%
    select(Sample1, Sample2, Matches, Mismatches, Jaccard.Index, Mash.like.distance, SNPs, SNP.distance)
 }

# colorvector=c("#009292","#b66dff","#db6d00","#6db6ff","#924900","#20db20","#ff6db6","#490092")
# colorvector <- c("#FFDAB9", "#FFA07A", "#FFC300", "#ef8f01", "#D2691E", "#8B4513", "#1C3A3A", "black")
colorvector <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#666666")

ska_matrix <- list()
ska_plot <- list()
for (dataset in names(ska_distance)) {

  other_half <- ska_distance[[dataset]] %>%
    select(Sample1, Sample2, SNP.distance) %>%
    mutate(temp_sample = Sample1) %>%
    mutate(Sample1 = Sample2,
           Sample2 = temp_sample) %>%
    select(-temp_sample)
  
  ska_matrix[[dataset]] <- ska_distance[[dataset]] %>%
    select(Sample1, Sample2, SNP.distance) %>%
    bind_rows(other_half) %>%
    pivot_wider(values_from = SNP.distance, names_from = Sample1, values_fill = 0) %>%
    arrange(Sample2) %>%
    column_to_rownames("Sample2")
  
  pcoa <- ska_matrix[[dataset]] %>%
    as.dist() %>%
    pcoa()
  
  perc_var_axis1 <- round(pcoa$values$Rel_corr_eig[1]*100, digits=1)
  perc_var_axis2 <- round(pcoa$values$Rel_corr_eig[2]*100, digits=1)
  
  plot_df <- pcoa$vectors %>%
    as.data.frame() %>%
    select(Axis.1, Axis.2) %>%
    rownames_to_column("Bee_pool") %>% 
    left_join(., metadata, by = "Bee_pool") %>%
    # select(Bee_pool, Axis.1, Axis.2, Country, Season, Health) %>%
    select(Bee_pool, Axis.1, Axis.2, Country, Season, Health, Mapped_reads) %>%
    group_by(Bee_pool) %>%
    mutate(Mapped_reads = sum(Mapped_reads)) %>%
    ungroup() %>%
    distinct() %>%
    mutate(Country = factor(Country, levels = c("BE", "CH", "DE", "FR", "NL", "PT", "RO", "UK")))
  
  ska_plot[[dataset]] <- ggplot(plot_df, aes(x=Axis.1, y=Axis.2, color=Country)) +
    labs(x = paste0("Axis.1 [", perc_var_axis1, "%]"), y=paste0("Axis.2 [", perc_var_axis2, "%]"), 
         title= paste0("SKA2 - ", dataset)) +
    geom_point() +
    scale_color_manual(values = colorvector) +
    theme_bw()
}

pcoa_wrap <- wrap_plots(ska_plot) +
  plot_layout(guides = "collect")
pcoa_wrap

ska_plot$phages +
  aes(color=Mapped_reads) +
  scale_color_continuous()

clustered <- list()
for (dataset in names(ska_matrix)) {
  clustered[[dataset]] <- ska_matrix[[dataset]] %>%
    as.dist() %>%
    hclust(method = "ward.D") %>%
    as.phylo()
}

# SNP mismatch trees
tree_plots <- list()
for (dataset in names(clustered)) {
  # tree_plots[[dataset]] <- ggtree(clustered[[dataset]], branch.length="none") +
  tree_plots[[dataset]] <- ggtree(clustered[[dataset]]) +
    # geom_tiplab(aes(color = substr(label, 1, 2))) +
    geom_tippoint(aes(color = substr(label, 1, 2))) +
    scale_color_manual(values = colorvector, name = "Country") +
    guides(color = guide_legend(override.aes = list(label = "•", size = 5))) +
    ggtitle(dataset)
}

tree_wrap <- wrap_plots(tree_plots) +
  plot_layout(guides = "collect")
tree_wrap

#########

pcoa_wrap <- wrap_plots(ska_plot$host / ska_plot$phages / ska_plot$MAG_mapped) +
  plot_layout(guides = "collect")
ggsave("output/SNP_analysis/SNP_PCoAs.pdf", tree_wrap,
       width = 7, height = 15)

tree_wrap <- wrap_plots(tree_plots$host / tree_plots$phages / tree_plots$MAG_mapped) +
  plot_layout(guides = "collect")
ggsave("output/SNP_analysis/SNP_trees.pdf", pcoa_wrap,
       width = 7, height = 15)

### Tests

TreeDistance(clustered$host, clustered$phages)
TreeDistance(clustered$host, clustered$MAG_mapped)
TreeDistance(clustered$phages, clustered$MAG_mapped)

mantel.test(ska_matrix$host, ska_matrix$phages)
mantel.test(ska_matrix$host, ska_matrix$MAG_mapped)
mantel.test(ska_matrix$phages, ska_matrix$MAG_mapped)

combinations <- combn(names(ska_matrix), 2)
mantels <- list()
for (i in 1:ncol(combinations)) {
  element1 <- combinations[1, i]
  element2 <- combinations[2, i]
  if (element1 %in% c("Acer_CO1", "Acer_genome", "Acer_mt", "host") & element2 %in% c("Acer_CO1", "Acer_genome", "Acer_mt", "host")) {
    next
  }
  if (element1 %in% c("all_terLs", "core", "core_terLs", "phages") & element2 %in% c("all_terLs", "core", "core_terLs", "phages")) {
    next
  }
  common_pools <- intersect(rownames(ska_matrix[[element1]]), rownames(ska_matrix[[element2]]))
  
  m1_filtered <- ska_matrix[[element1]][common_pools, common_pools, drop = FALSE]
  m2_filtered <- ska_matrix[[element2]][common_pools, common_pools, drop = FALSE]
  
  mantels[[paste0(element1, "_", element2)]] <- mantel.test(m1_filtered, m2_filtered)
}

for (test in names(mantels)) {
  if (paste0(mantels[[test]]$p <= 0.05)) {
    print(paste0(mantels[[test]]$p, " - ", test))
  }
}





#
