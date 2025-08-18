library(ape)
library(patchwork)
library(ggtree)
library(RLdbRDA)
library(vegan)
library(tidyverse)

set.seed(1)

source("scripts/R/custom_rldbrda.R")
source("scripts/R/helpers/mixed_helpers.R")


metadata <- readRDS("data/metadata.RDS")

relevant_meta <- metadata %>%
  select(Bee_pool, Country) %>%
  filter(!is.na(Country)) %>%
  distinct()

datasets <- c("bee", "phages", "bacteria")

raw_input <- list()

for (dataset in c("bee", "phages", "bacteria")) {
  raw_input[[dataset]] <- read.delim(paste0("output/SNP_analysis/SKA_SNP_distances/", dataset, ".distances.tsv")) 
}

ska_distance <- list()
for (dataset in datasets) {
  ska_distance[[dataset]] <- raw_input[[dataset]] %>% tibble() %>%
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

colorvector <- list()
colorvector$phages <- country_colors
colorvector$bee <- country_colors
colorvector$bacteria <- country_colors

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
    select(Bee_pool, Axis.1, Axis.2, Country) %>%
    distinct()
  
  centroids <- plot_df %>%
    group_by(Country) %>%
    summarise(centroid.1 = mean(Axis.1), centroid.2 = mean(Axis.2))
    
  ska_plot[[dataset]] <- ggplot(plot_df, aes(x = Axis.1, y = Axis.2, color = Country)) +
    labs(x = paste0("Axis.1 [", perc_var_axis1, "%]"), y=paste0("Axis.2 [", perc_var_axis2, "%]"), 
         title = paste0("SKA - ", dataset)) +
    geom_point() +
    
    # Centroids and circles around groups (clutters the plots)
    # geom_point(data = centroids, aes(x = centroid.1, y = centroid.2), shape=13, size=4) +
    # stat_ellipse(
    #   # aes(group = Variable1), # Group by Variable1
    #   aes(group = Country), # Group by Country
    #   level = 0.95, type = "t", linetype = "dashed"
    # ) +
    
    scale_color_manual(values = colorvector[[dataset]]) +
    theme_bw()
  ska_plot[[dataset]]
}

pcoa_wrap <- wrap_plots(ska_plot) +
  plot_layout(guides = "collect")

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
    # geom_tiplab(aes(label = Country), hjust = -0.3, size = 2) +
    geom_tippoint(aes(color = Country)) +
    scale_color_manual(values = colorvector[[dataset]], name = "Country") +
    guides(color = guide_legend(override.aes = list(label = "•", size = 5))) +
    ggtitle(dataset) +
    coord_cartesian(clip = "off")
}

tree_wrap <- wrap_plots(tree_plots) +
  plot_layout(guides = "collect")

#########
# RDAs

custom_colors <- c("#8B4513", "#1C3A3A", "#FFA07A", "#1C3A3A")

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
  
  p <- RLdbRDA::plot_dbrda(plot_data)
  p$data$bar_color <- custom_colors
  
  RDA_plots[[set]] <- p + aes(fill = bar_color) + scale_fill_identity()
}

common_legend <- legend_factory(title = "Thing", 
                                items = c("R2_A", "R2_B", "Cumulative R2"),
                                colors = c("#8B4513", "#FFA07A", "#1C3A3A"),
                                position = "bottom")

RDA_patch_horizontal <- (RDA_plots$bee + RDA_plots$bacteria + RDA_plots$phages + plot_layout(axes = "collect")) /
  common_legend +
  plot_layout(heights = c(12,1))

RDA_patch_vertical <- RDA_plots$bee / RDA_plots$bacteria / RDA_plots$phages / common_legend +
  plot_layout(axes = "collect", heights = c(rep(10,3), 1))

#####
### Tests

mantel_tests <- list()
for (mant in c("bee_phages", "bee_bacteria", "phages_bacteria")) {
  first_set <- str_split(mant, "_")[[1]][[1]]
  second_set <- str_split(mant, "_")[[1]][[2]]
  
  mantel_tests[[mant]] <- mantel(ska_matrix[[first_set]], ska_matrix[[second_set]], method = "spearman") %>%
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

#####
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

ggsave("output/R/SNP_analysis/SNP_RDA_horizontal.pdf", RDA_patch_horizontal,
       width = 11, height = 3)
ggsave("output/R/SNP_analysis/SNP_RDA_vertical.pdf", RDA_patch_vertical,
       width = 8, height = 4)

for (dataset in names(ska_plot)) {
  ggsave(paste0("output/R/SNP_analysis/SNP_PCoA_", dataset, ".pdf"),
         ska_plot[[dataset]], width = 6, height = 5)
  ggsave(paste0("output/R/SNP_analysis/SNP_tree_", dataset, ".pdf"),
         tree_plots[[dataset]], width = 7, height = 10)
  ggsave(paste0("output/R/SNP_analysis/SNP_RDA_", dataset, ".pdf"),
         RDA_plots[[dataset]], width = 5, height = 4)
  write_delim(RDAs[[dataset]], 
              paste0("output/R/SNP_analysis/SNP_RDA_", dataset, ".tsv"),
              delim = "\t")
}

for (test in names(mantel_tests)) {
  write_delim(mantel_tests[[test]],
              paste0("output/R/SNP_analysis/SNP_mantel_", test, ".tsv"),
              delim = "\t", col_names = TRUE)
}
