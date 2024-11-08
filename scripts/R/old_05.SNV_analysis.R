library(tidyverse)
library(ape)

concatenated_core_terLs.filtered.mann <- read.delim("output/SNV_analysis/metaSNV_terL_core/distances/concatenated_core_terLs.filtered.mann.dist") %>%
  column_to_rownames("X") %>%
  as.dist()

metadata <- readRDS("output/R/R_variables/metadata.RDS")




concenated_terLs_pcoa <- pcoa(concatenated_core_terLs.filtered.mann)

perc_var_axis1 <- round(concenated_terLs_pcoa$values$Relative_eig[1]*100, digits=1)
perc_var_axis2 <- round(concenated_terLs_pcoa$values$Relative_eig[2]*100, digits=1)
# perc_var_axis3 <- round(concenated_terLs_pcoa$values$Relative_eig[3]*100, digits=1)
# perc_var_axis4 <- round(concenated_terLs_pcoa$values$Relative_eig[4]*100, digits=1)


colorvector=c("#009292","#b66dff","#db6d00","#6db6ff","#924900","#20db20","#ff6db6","#490092")

plot_df <- concenated_terLs_pcoa$vectors %>%
  as.data.frame() %>%
  select(Axis.1, Axis.2) %>%
  rownames_to_column("bam") %>% 
  mutate(Bee_pool = sub("^(([^_]*_){2}[^_]*).*", "\\1", bam)) %>% 
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

ggplot(plot_df, aes(x=Axis.1, y=Axis.2, color=Country)) +
  labs(x = paste0("Axis.1 [", perc_var_axis1, "%]"), y=paste0("Axis.2 [", perc_var_axis2, "%]"), 
       title= "core terL Manhattan") +
  geom_point() +
  scale_color_manual(values = colorvector) +
  theme_bw()
