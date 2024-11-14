library(ape)
library(patchwork)
library(ggtree)
library(TreeDist)
library(tidyverse)

set.seed(1)

metadata <- readRDS("output/R/R_variables/metadata.RDS")
# datasets <- c("Acer_CO1", "Acer_genome", "Acer_mt", "host", "phages", "MAG_mapped", "core", "all_terLs", "core_terLs")
datasets <- c("Acer_CO1", "all_terLs", "host", "phages")

ska2_distance <- list()
ska2_tree <- list()
ska2_raxml_besttree <- list()
ska2_raxml_resulttree <- list()
for (dataset in datasets) {
  ska2_distance[[dataset]] <- read.delim(paste0("output/SNP_analysis/ska2_", dataset, "/ska2_", dataset,".distance"))
  # ska2_tree[[dataset]] <- read.tree(paste0("output/SNP_analysis/ska2_", dataset, "/ska2_", dataset,".tree"))
  # ska2_raxml_besttree[[dataset]] <- read.tree(paste0("output/SNP_analysis/ska2_", dataset, "/RAxML_bestTree.ska2_", dataset,".raxml.tree"))
  # ska2_raxml_resulttree[[dataset]] <- read.tree(paste0("output/SNP_analysis/ska2_", dataset, "/RAxML_result.ska2_", dataset,".raxml.tree"))
}

# colorvector=c("#009292","#b66dff","#db6d00","#6db6ff","#924900","#20db20","#ff6db6","#490092")
# colorvector <- c("#FFDAB9", "#FFA07A", "#FFC300", "#ef8f01", "#D2691E", "#8B4513", "#1C3A3A", "black")
colorvector <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#666666")

ska2_matrix <- list()
ska2_plot <- list()
snp_matrix <- list()
for (dataset in names(ska2_distance)) {

  other_half <- ska2_distance[[dataset]] %>%
    select(-Distance) %>%
    mutate(temp_sample = Sample1) %>%
    mutate(Sample1 = Sample2,
           Sample2 = temp_sample) %>%
    select(-temp_sample)
  
  ska2_matrix[[dataset]] <- ska2_distance[[dataset]] %>%
    select(-Distance) %>%
    bind_rows(other_half) %>%
    pivot_wider(values_from = Mismatches, names_from = Sample2, values_fill = 0) %>%
    column_to_rownames("Sample1")
  ska2_matrix[[dataset]] <- ska2_matrix[[dataset]][match(colnames(ska2_matrix[[dataset]]), rownames(ska2_matrix[[dataset]])), ]
  
  pcoa <- ska2_matrix[[dataset]] %>%
    as.dist() %>%
    pcoa()
  
  perc_var_axis1 <- round(pcoa$values$Rel_corr_eig[1]*100, digits=1)
  perc_var_axis2 <- round(pcoa$values$Rel_corr_eig[2]*100, digits=1)
  
  plot_df <- pcoa$vectors %>%
    as.data.frame() %>%
    select(Axis.1, Axis.2) %>%
    rownames_to_column("Bee_pool") %>% 
    left_join(., metadata, by = "Bee_pool") %>%
    select(Bee_pool, Axis.1, Axis.2, Country, Season, Health) %>%
    distinct() %>%
    mutate(Country = factor(Country, levels = c("BE", "CH", "DE", "FR", "NL", "PT", "RO", "UK")))
  
  ska2_plot[[dataset]] <- ggplot(plot_df, aes(x=Axis.1, y=Axis.2, color=Country)) +
    labs(x = paste0("Axis.1 [", perc_var_axis1, "%]"), y=paste0("Axis.2 [", perc_var_axis2, "%]"), 
         title= paste0("SKA2 - ", dataset)) +
    geom_point() +
    scale_color_manual(values = colorvector) +
    theme_bw()
}

wrap <- wrap_plots(ska2_plot) +
  plot_layout(guides = "collect")
wrap

combinations <- combn(names(ska2_matrix), 2)
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
  common_pools <- intersect(rownames(ska2_matrix[[element1]]), rownames(ska2_matrix[[element2]]))

  m1_filtered <- ska2_matrix[[element1]][common_pools, common_pools, drop = FALSE]
  m2_filtered <- ska2_matrix[[element2]][common_pools, common_pools, drop = FALSE]

  mantels[[paste0(element1, "_", element2)]] <- mantel.test(m1_filtered, m2_filtered)
}

for (test in names(mantels)) {
  if (paste0(mantels[[test]]$p <= 0.05)) {
    print(paste0(mantels[[test]]$p, " - ", test))
  }
}


clustered <- list()
for (dataset in names(ska2_matrix)) {
  clustered[[dataset]] <- ska2_matrix[[dataset]] %>%
    as.dist() %>%
    hclust(method = "ward.D") %>%
    as.phylo()
}

# SNP mismatch trees
tree_plots <- list()
for (dataset in names(clustered)) {
  # tree_plots[[dataset]] <- ggtree(clustered[[dataset]], branch.length="none") +
  tree_plots[[dataset]] <- ggtree(clustered[[dataset]]) +
      geom_tiplab(aes(color = substr(label, 1, 2))) +
    scale_color_manual(values = colorvector, name = "Country") +
    guides(color = guide_legend(override.aes = list(label = "•", size = 10))) +
    ggtitle(dataset)
}
tree_plots$phages
tree_plots$host
tree_plots$Acer_genome
tree_plots$Acer_CO1
tree_plots$all_terLs
TreeDistance(clustered$Acer_CO1, clustered$all_terLs)

# Alignment trees
for (dataset in names(ska2_tree)) {
  ggtree(ska2_tree[[dataset]]) +
    geom_tiplab(aes(color = substr(label, 1, 2))) +
    scale_color_manual(values = colorvector, name = "Country") +
    guides(color = guide_legend(override.aes = list(label = "•", size = 10))) +
    ggtitle(dataset)
  
  ggtree(ska2_raxml_besttree[[dataset]]) +
    geom_tiplab(aes(color = substr(label, 1, 2))) +
    scale_color_manual(values = colorvector, name = "Country") +
    guides(color = guide_legend(override.aes = list(label = "•", size = 10))) +
    ggtitle(dataset)
  
  ggtree(ska2_raxml_resulttree[[dataset]]) +
    geom_tiplab(aes(color = substr(label, 1, 2))) +
    scale_color_manual(values = colorvector, name = "Country") +
    guides(color = guide_legend(override.aes = list(label = "•", size = 10))) +
    ggtitle(dataset)
}


# system("mkdir -p output/R/SNP_analysis")
# ggsave("output/R/SNP_analysis/SKA2.pdf", wrap,
#        width = 16, height = 15)



# q1 <- matrix(runif(22500), nrow = 150)
# q2 <- matrix(runif(22500), nrow = 150)
# diag(q1) <- diag(q2) <- 0
# mantel.test(q1, q2, graph = TRUE, nperm = 10000,
#             main = "Mantel test: a random example with 150 X 150 matrices
# representing asymmetric relationships",
#             xlab = "z-statistic", ylab = "Density",
#             sub = "The vertical line shows the observed z-statistic")


bla <- c("RO_17381_sum", "RO_17380_sum", "BE_16562_spr", "BE_16562_aut")
ska2_distance$Acer_CO1 %>%
  tibble() %>%
  filter(Sample1 %in% bla) %>%
  filter(Sample2 %in% bla) %>%
  select(-Distance) %>%
  pivot_wider(names_from = Sample2, values_from = Mismatches)







#
