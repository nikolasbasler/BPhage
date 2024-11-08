library(tidyverse)
library(ape)
library(patchwork)

metadata <- readRDS("output/R/R_variables/metadata.RDS")

ska.distances <- list()
ska.distances$terLs <- read.delim("output/SNP_analysis/from_Lander/ska_plots/terLs/ska.distances.tsv")
ska.distances$host <- read.delim("output/SNP_analysis/from_Lander/ska_plots/host/ska.distances.tsv")

ska1_matrix <- list()
ska1_plot <- list()
for (dataset in names(ska.distances)) {
  other_half <- ska.distances[[dataset]] %>%
    select(Sample.1, Sample.2, SNP.distance) %>%
    mutate(temp_sample = Sample.1) %>%
    mutate(Sample.1 = Sample.2,
           Sample.2 = temp_sample) %>%
    select(-temp_sample)
  
  
  ska1_matrix[[dataset]] <- ska.distances[[dataset]] %>%
    select(Sample.1, Sample.2, SNP.distance) %>%
    bind_rows(other_half) %>%
    pivot_wider(values_from = SNP.distance, names_from = Sample.2, values_fill = 0) %>%
    column_to_rownames("Sample.1")
  
  ska_pcoa <- ska1_matrix[[dataset]] %>%
    as.dist() %>%
    pcoa()
  
  perc_var_axis1 <- round(ska_pcoa$values$Rel_corr_eig[1]*100, digits=1)
  perc_var_axis2 <- round(ska_pcoa$values$Rel_corr_eig[2]*100, digits=1)
  
  colorvector=c("#009292","#b66dff","#db6d00","#6db6ff","#924900","#20db20","#ff6db6","#490092")
  
  plot_df <- ska_pcoa$vectors %>%
    as.data.frame() %>%
    select(Axis.1, Axis.2) %>%
    rownames_to_column("Bee_pool") %>% 
    left_join(., metadata, by = "Bee_pool") %>%
    select(Bee_pool, Axis.1, Axis.2, Country, Season, Health) %>%
    distinct() %>%
    mutate(Country = factor(Country, levels = c("BE", "CH", "DE", "FR", "NL", "PT", "RO", "UK")))
  
  
  plot_df %>%
    group_by(Country) %>%
    count()
  plot_df %>%
    group_by(Season) %>%
    count()
  plot_df %>%
    group_by(Health) %>%
    count()
  
  ska1_plot[[dataset]] <- ggplot(plot_df, aes(x=Axis.1, y=Axis.2, color=Country)) +
    labs(x = paste0("Axis.1 [", perc_var_axis1, "%]"), y=paste0("Axis.2 [", perc_var_axis2, "%]"), 
         title = paste0("SNP distance SKA1 - ", dataset)) +
    geom_point() +
    scale_color_manual(values = colorvector) +
    theme_bw()
}

mantel.test(ska1_matrix$terLs, ska1_matrix$host, graph = TRUE)

# wrap <- wrap_plots(ska1_plot) +
#   plot_layout(guides = "collect")
# system("mkdir -p output/R/SNP_analysis")
# ggsave("output/R/SNP_analysis/SKA1.pdf", wrap,
#        width = 12, height = 5.5)
