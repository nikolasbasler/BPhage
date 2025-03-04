library(ape)
library(patchwork)
library(ggtree)
library(RLdbRDA)
source("scripts/R/custom_rldbrda.R")
library(vegan)
library(tidyverse)

set.seed(1)

metadata <- readRDS("output/R/R_variables/metadata.RDS")
SKA_subspecies_SRAs <- read.delim("data/SKA_subspecies_SRAs.tsv") %>%
  select(-Pool_name) %>%
  mutate(Country = NA) %>%
  rename(Bee_pool = Accession)

relevant_meta <- metadata %>%
  tibble() %>%
  select(Bee_pool, Country) %>%
  mutate(Lineage = NA, 
         Subspecies = NA) %>%
  rbind(SKA_subspecies_SRAs) %>%
  filter(!is.na(Bee_pool)) %>%
  mutate(Lineage = factor(Lineage, levels = c("A", "M", "O", "C")),
         Subspecies = factor(Subspecies, levels = c("ruttneri", "iberiensis", "mellifera", 
                                                    "anatoliaca", "caucasia", "remipes", 
                                                    "cypria", "ligustica", "carnica", "adami", 
                                                    "cecropia", "macedonica", "carpatica", 
                                                    "rodopica"))) %>%
  mutate(Variable1 = ifelse(is.na(Country), as.character(Subspecies), as.character(Country))) %>%
  mutate(Variable1 = factor(Variable1, levels = union(levels(Country), levels(Subspecies))),
         Variable1 = droplevels(Variable1)) %>%
  mutate(Variable2 = ifelse(is.na(Country), as.character(Lineage), as.character(Country))) %>%
  mutate(Variable2 = factor(Variable2, levels = union(levels(Country), levels(Lineage))),
         Variable2 = droplevels(Variable2))

# datasets <- c("bee", "subspecies", "bee_and_subspecies", "phages", "bacteria")
datasets <- c("bee", "phages", "bacteria")

raw_input <- list()
# for (dataset in c("bee", "subspecies", "bee_and_subspecies")) {
#   raw_input[[dataset]] <- read.delim(paste0("output/SNP_analysis/SKA_SNP_distances_region_mapped/", dataset, ".distances.tsv")) %>%
#     filter(!str_starts(Sample.1, "FR_19771_aut")) %>%
#     filter(!str_starts(Sample.2, "FR_19771_aut")) %>%
#     filter(!str_starts(Sample.1, "BE_16558_sum")) %>%
#     filter(!str_starts(Sample.2, "BE_16558_sum")) %>%
#     filter(!str_starts(Sample.1, "PT_19409_sum")) %>%
#     filter(!str_starts(Sample.2, "PT_19409_sum")) %>%
#     filter(!str_starts(Sample.1, "PT_19414_aut")) %>%
#     filter(!str_starts(Sample.2, "PT_19414_aut")) %>%
#     filter(!str_starts(Sample.1, "PT_19409_spr")) %>%
#     filter(!str_starts(Sample.2, "PT_19409_spr")) # Removing a few weird outliers. The mantel tests would be the same with them.
# }
for (dataset in c("bee", "phages", "bacteria")) {
  raw_input[[dataset]] <- read.delim(paste0("output/SNP_analysis/SKA_SNP_distances/", dataset, ".distances.tsv")) 
}

ska_distance <- list()
for (dataset in datasets) {
  ska_distance[[dataset]] <- raw_input[[dataset]] %>%
    mutate(Sample.1 = str_replace(Sample.1, "subspecies", "NA_NA_subpsecies")) %>%
    mutate(Sample.2 = str_replace(Sample.2, "subspecies", "NA_NA_subpsecies")) %>%
    separate_wider_delim(Sample.1, "_", names = c("country", "hive", "season"), too_many = "drop") %>%
    mutate(Sample1 = paste(country, hive, season, sep="_")) %>%
    select(-country, -hive, -season) %>%
    separate_wider_delim(Sample.2, "_", names = c("country", "hive", "season"), too_many = "drop") %>%
    mutate(Sample2 = paste(country, hive, season, sep="_")) %>%
    select(Sample1, Sample2, Matches, Mismatches, Jaccard.Index, Mash.like.distance, SNPs, SNP.distance) %>%
    mutate(Sample1 = str_replace(Sample1, "_NA_NA", "")) %>%
    mutate(Sample2 = str_replace(Sample2, "_NA_NA", "")) %>% 
    filter(!is.na(SNP.distance))
}

country_colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#666666")
# subspec_colors <- c("#FF0000", "#00FF00", "#0000FF", "#FF00FF", "#00FFFF", "#FFFF00", "#800000", "#808000", "#008080", "#800080", "#FFA500", "#A52A2A", "#8B4513", "#4B0082")
# lineage_colors <- c("#FF0000", "#00FF00", "#0000FF", "#FF00FF")
# lineage_colors2 <- c("#FFA500", "#00FFFF", "#FFFF00", "#BFFF00")

colorvector <- list()
colorvector$phages <- country_colors
colorvector$bee <- country_colors
colorvector$bacteria <- country_colors
# colorvector$subspecies <- subspec_colors
# colorvector$subspecies_lin <- lineage_colors
# colorvector$bee_and_subspecies <- c(country_colors, subspec_colors)
# colorvector$bee_and_subspecies_lin <- c(country_colors, lineage_colors)

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
    left_join(., relevant_meta, by = "Bee_pool") %>%
    select(Bee_pool, Axis.1, Axis.2, Variable1, Variable2) %>%
    distinct()
  
  centroids <- plot_df %>%
    group_by(Variable1) %>%
    summarise(centroid.1 = mean(Axis.1), centroid.2 = mean(Axis.2))
    
  ska_plot[[dataset]] <- ggplot(plot_df, aes(x = Axis.1, y = Axis.2, color = Variable1)) +
    labs(x = paste0("Axis.1 [", perc_var_axis1, "%]"), y=paste0("Axis.2 [", perc_var_axis2, "%]"), 
         title = paste0("SKA - ", dataset)) +
    geom_point() +
    
    # Centroids and circles around groups
    # geom_point(data = centroids, aes(x = centroid.1, y = centroid.2), shape=13, size=4) +
    # stat_ellipse(
    #   aes(group = Variable1), # Group by Variable1
    #   level = 0.95, type = "t", linetype = "dashed"
    # ) +
    
    scale_color_manual(values = colorvector[[dataset]]) +
    theme_bw()
  ska_plot[[dataset]]

  if (dataset %in% c("subspecies", "bee_and_subspecies")) {
    ska_plot[[paste0(dataset, "_lineage")]] <- ska_plot[[dataset]] +
      aes(color = Variable2) +
      scale_color_manual(values = colorvector[[paste0(dataset, "_lin")]]) +
      labs(title = paste0("SKA - ", dataset, "_lineage"))
  }
}

pcoa_wrap <- wrap_plots(ska_plot) +
  plot_layout(guides = "collect")
pcoa_wrap

#### Trees

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
  tree_plots[[dataset]] <- ggtree(clustered[[dataset]]) %<+% relevant_meta +
    geom_tiplab(aes(label = Variable1), hjust = -0.3, size = 2) +
    geom_tippoint(aes(color = Variable1)) +
    scale_color_manual(values = colorvector[[dataset]], name = "Variable1") +
    guides(color = guide_legend(override.aes = list(label = "•", size = 5))) +
    ggtitle(dataset) +
    coord_cartesian(clip = "off")
  # if (dataset %in% c("subspecies", "bee_and_subspecies")) {
  if (dataset == "subspecies") {
      tree_plots[[dataset]] <- tree_plots[[dataset]] +
      geom_tiplab(aes(label = Variable1, color = Variable2), hjust = -0.3, size = 2) +
      scale_color_manual(values =
                           c(colorvector[[dataset]], lineage_colors2))
  }
  
  if (dataset == "bee_and_subspecies") {
    tree_plots[[dataset]] <- tree_plots[[dataset]] +
      geom_tiplab(aes(label = Variable1, color = Variable2), hjust = -0.3, size = 2) +
      scale_color_manual(values =
                           c(rep("black", 8), subspec_colors, lineage_colors2))
  }
}

tree_wrap <- wrap_plots(tree_plots) +
  plot_layout(guides = "collect")
tree_wrap

#########
# RDAs

RDAs <- list()
RDA_plots <- list()
for (set in names(ska_matrix)) {
  dist_matrix <- ska_matrix[[set]] %>% 
    as.dist()
  
  meta_filt <- metadata %>%
    filter(Bee_pool %in% rownames(ska_matrix[[set]])) %>%
    rownames_to_column("bla") %>%
    select(Bee_pool, Season, Country) %>%
    distinct() %>%
    column_to_rownames("Bee_pool")
 
  r_names_dist <- rownames(ska_matrix[[set]])
  if (!identical(r_names_dist, rownames(meta_filt))) {
    print("ROWNAMES OF DIST MATRIX AND METADATA DONT MATCH! ABORTING")
    break
  }
  RDAs[[set]] <- custom_rldbrda(dist_matrix, meta_filt)
  # detach(meta)
  
  plot_data <- RLdbRDA::prepare_plot_data(RDAs[[set]])
  
  RDA_plots[[set]] <- RLdbRDA::plot_dbrda(plot_data) + 
    scale_fill_manual(values=c("#ef8f01", "#8B4513"),
                      labels=c(bquote(R^2), bquote('Cumulative' ~ R^2))) +
    ggtitle(paste0(set, " SNP distance")) +
    theme(legend.position = "bottom")
}

RDA_patch_vertical <- RDA_plots$bee + RDA_plots$phages + RDA_plots$bacteria +
  plot_layout(guides = 'collect', axes = "collect")  & theme(legend.position = "bottom")
RDA_patch_horizontal <- RDA_plots$bee / RDA_plots$phages / RDA_plots$bacteria +
  plot_layout(guides = 'collect', axes = "collect")  & theme(legend.position = "bottom")

#####

### Tests
# 
# common_bee_pools <- intersect(colnames(ska_matrix$bee), colnames(ska_matrix$phages)) %>%
#   intersect(colnames(ska_matrix$bacteria))
# 
# ska_matrix_filt <- list()
# 
# ska_matrix_filt$bee <- ska_matrix$bee %>%
#   select(any_of(common_bee_pools)) %>%
#   rownames_to_column("pools") %>%
#   filter(pools %in% common_bee_pools) %>% 
#   column_to_rownames("pools")
# ska_matrix_filt$phages <- ska_matrix$phages %>%
#   select(any_of(common_bee_pools)) %>%
#   rownames_to_column("pools") %>%
#   filter(pools %in% common_bee_pools) %>% 
#   column_to_rownames("pools")
# ska_matrix_filt$bacteria <- ska_matrix$bacteria %>%
#   select(any_of(common_bee_pools)) %>%
#   rownames_to_column("pools") %>%
#   filter(pools %in% common_bee_pools) %>% 
#   column_to_rownames("pools")

#####
# TreeDistance(clustered$bee, clustered$phages)
# TreeDistance(clustered$bee, clustered$bacteria)
# TreeDistance(clustered$phages, clustered$bacteria)

mantel_tests <- list()
for (mant in c("bee_phages", "bee_bacteria", "phages_bacteria")) {
  first_set <- str_split(mant, "_")[[1]][[1]]
  second_set <- str_split(mant, "_")[[1]][[2]]
  
  mantel_tests[[mant]] <- mantel(ska_matrix[[first_set]], ska_matrix[[second_set]]) %>%
    unlist() %>%
    t() %>% 
    as.data.frame() %>%
    mutate(coeff = as.numeric(statistic),
           p = as.numeric(signif),
           method = unlist(method)) %>%
    select(method, coeff, p)
}


pcoa_wrap_vertical <- wrap_plots(ska_plot$bee / ska_plot$phages / ska_plot$bacteria) +
  plot_layout(guides = "collect")
pcoa_wrap_horizontal <- wrap_plots(ska_plot$bee + ska_plot$phages + ska_plot$bacteria) +
  plot_layout(guides = "collect")
tree_wrap_vertical <- wrap_plots(tree_plots$bee / tree_plots$phages / tree_plots$bacteria) +
  plot_layout(guides = "collect")
tree_wrap_horizontal <- wrap_plots(tree_plots$bee + tree_plots$phages + tree_plots$bacteria) +
  plot_layout(guides = "collect")

# Save files

system("mkdir -p output/R/SNP_analysis/")
ggsave("output/R/SNP_analysis/SNP_PCoAs_vertical.pdf", pcoa_wrap_vertical,
       width = 7, height = 15)
ggsave("output/R/SNP_analysis/SNP_PCoAs_horizontal.pdf", pcoa_wrap_horizontal,
       width = 16, height = 5)
ggsave("output/R/SNP_analysis/SNP_trees_vertical.pdf", tree_wrap_vertical,
       width = 7, height = 20)
ggsave("output/R/SNP_analysis/SNP_trees_horizontal.pdf", tree_wrap_horizontal,
       width = 15, height = 10)

ggsave("output/R/SNP_analysis/SNP_RDA_vertical.pdf", RDA_patch_vertical,
       width = 9, height = 3)
ggsave("output/R/SNP_analysis/SNP_RDA_horizontal.pdf", RDA_patch_horizontal,
       width = 8, height = 4)



for (dataset in names(ska_plot)) {
  ggsave(paste0("output/R/SNP_analysis/SNP_PCoA_", dataset, ".pdf"),
         ska_plot[[dataset]], width = 6, height = 5)
  ggsave(paste0("output/R/SNP_analysis/SNP_tree_", dataset, ".pdf"),
         tree_plots[[dataset]], width = 7, height = 10)
  ggsave(paste0("output/R/SNP_analysis/SNP_RDA_", dataset, ".pdf"),
         RDA_plots[[dataset]], width = 5, height = 4)
}


for (test in names(mantel_tests)) {
  write_delim(mantel_tests[[test]],
              paste0("output/R/SNP_analysis/SNP_mantel_", test, ".tsv"),
              delim = "\t", col_names = TRUE)
}
