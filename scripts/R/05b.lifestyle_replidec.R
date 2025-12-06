library(tidyverse)

classification <- readRDS("data/classification.RDS")
lifestyle_colors <- c("Virulent" = "#1C3A3A", "Chronic" =  "#8B4513", "Temperate" ="#FFC300")

present_in_all_countries <- read_lines("data/core_contigs.txt")

BC_predict <- read.delim("output/lifestyle/replidec/BC_predict.summary") %>%
  tibble()

## PLOT

updated_classif <- classification %>%
  left_join(., BC_predict[c("sample_name", "final_label")], by = join_by("contig" == "sample_name")) %>%
  mutate(Lifestyle_replidec = final_label) %>% 
  relocate(Lifestyle_replidec, .after = Core) %>%
  select(-final_label)

replidec_lifestyle_tibble <- updated_classif %>%
  filter(
    Class == "Caudoviricetes",
    # Order == "Microviruses",
  ) %>%
  count(Core, Lifestyle_replidec) %>%
  group_by(Core) %>%
  mutate(
    prop = n / sum(n),
    Core = factor(Core, levels = c("yes", "no")),
    Lifestyle_replidec = factor(Lifestyle_replidec, levels = c("Virulent", "Chronic", "Temperate"))
  ) %>%
  ungroup()

replidec_lifestyle_plot <- replidec_lifestyle_tibble %>%
  ggplot(aes(y = Core, x = prop, fill = Lifestyle_replidec)) +
  geom_col(color = "black") +
  scale_fill_manual(values = lifestyle_colors) +
  theme_minimal() +
  theme(
    legend.position = "top", legend.justification = "center",
    legend.title = element_blank()
    ) +
  guides(fill = guide_legend(reverse=TRUE))

##### Save files

system("mkdir -p output/R/lifestyle")
ggsave("output/R/lifestyle/replidec_horizontal_all_caudos.pdf",
       replidec_lifestyle_plot, 
       width = 5.25, height = 2.81
)
write_delim(replidec_lifestyle_tibble, "output/R/lifestyle/replidec_horizontal_all_caudos.tsv",
            delim = "\t")

## For convenience and to avoid back-tracking:
# write_csv(updated_classif, "output/R/classification.csv")

