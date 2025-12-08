library(tidyverse)

classification <- readRDS("data/classification.RDS")

lifestyle_colors <- c("Virulent" = "#1C3A3A", 
                      "uncertain_Virulent" = "#338080",
                      "Chronic" =  "#8B4513", 
                      "Temperate" ="#FFC300", 
                      "uncertain_Temperate" = "#FFDAB9")

bphage_and_extended.fasta <- read.delim("output/lifestyle/bacphlip/bphage_and_extended.fasta.bacphlip") %>%
  tibble() %>%
  rename(contig = X)

present_in_all_countries <- read_lines("data/core_contigs.txt")

threshold = 0.90

bacphlip_assignment <- classification %>%
  left_join(., bphage_and_extended.fasta, by = "contig") %>%
  select(contig, Virulent, Temperate) %>%
  mutate(assignment = case_when(
    Virulent >= threshold ~ "Virulent",
    Temperate >= threshold ~ "Temperate",
    Virulent > Temperate & Virulent < threshold ~ "uncertain_Virulent",
    Temperate > Virulent & Temperate < threshold ~ "uncertain_Temperate",
    .default = "uncertain"
  ))

## PLOT

updated_classif <- classification %>%
  left_join(., bacphlip_assignment, by = "contig") %>%
  mutate(Lifestyle_bacphlip = assignment) %>% 
  relocate(Lifestyle_bacphlip, .after = Core) %>%
  select(-c(assignment, Virulent, Temperate))

bacphlip_lifestyle_tibble <- updated_classif %>%
  filter(
    Class == "Caudoviricetes",
    # Order == "Microviruses"
  ) %>%
  count(Core, Lifestyle_bacphlip) %>%
  group_by(Core) %>%
  mutate(
    prop = n / sum(n),
    Core = factor(Core, levels = c("yes", "no")),
    Lifestyle_bacphlip = factor(Lifestyle_bacphlip, levels = c("Virulent", "uncertain_Virulent", "uncertain_Temperate", "Temperate"))
  ) %>%
  ungroup()

bacphlip_lifestyle_plot <- bacphlip_lifestyle_tibble %>%
  ggplot(aes(y = Core, x = prop, fill = Lifestyle_bacphlip)) +
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
ggsave("output/R/lifestyle/bacphlip_horizontal_all_caudos.pdf",
       bacphlip_lifestyle_plot, 
       width = 5.25, height = 2.81
)
write_delim(bacphlip_lifestyle_tibble, "output/R/lifestyle/bacphlip_horizontal_all_caudos.tsv",
            delim = "\t")

## For convenience and to avoid back-tracking:
# write_csv(updated_classif, "output/R/classification.csv")

