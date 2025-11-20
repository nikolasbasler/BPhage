library(tidyverse)

classification <- readRDS("data/classification.RDS")
lifestyle_colors <- c("Virulent" = "#1C3A3A", "Chronic" =  "#8B4513", "Temperate" ="#FFC300")

bphage_and_extended.fasta <- read.delim("output/lifestyle/bacphlip/bphage_and_extended.fasta.bacphlip") %>%
  tibble() %>%
  rename(contig = X)

present_in_all_countries <- read_lines("data/core_contigs.txt")


# 
# classification %>%
#   filter(
#     Class == "Caudoviricetes",
#     completeness == 100
#   ) %>%
#   count(Lifestyle_replidec, Lifestyle_bacphlip)
# 
# classification %>%
#   filter(
#     Class == "Caudoviricetes",
#     completeness == 100
#   ) %>%  group_by(Lifestyle_bacphlip) %>%
#   summarise(count = n())
# 
# repli <- list()
# repli$temp <- classification %>%
#   filter(
#     Class == "Caudoviricetes",
#     completeness == 100,
#     Lifestyle_replidec == "Temperate"
#     ) %>%
#   distinct(contig)
# repli$vir <- classification %>%
#   filter(
#     Class == "Caudoviricetes",
#     completeness == 100,
#     Lifestyle_replidec == "Virulent"
#   ) %>%
#   distinct(contig)
# 
# bacphlip <- list()
# bacphlip$temp <- classification %>%
#   filter(
#     Class == "Caudoviricetes",
#     completeness == 100,
#     Lifestyle_bacphlip == "Temperate"
#   ) %>%
#   distinct(contig)
# bacphlip$vir <- classification %>%
#   filter(
#     Class == "Caudoviricetes",
#     completeness == 100,
#     Lifestyle_bacphlip == "Virulent"
#   ) %>%
#   distinct(contig)

#

threshold = 0.90

bacphlip_assignment <- classification %>%
  left_join(., bphage_and_extended.fasta, by = "contig") %>%
  select(contig, Virulent, Temperate) %>%
  mutate(assignment = case_when(
    Virulent >= threshold ~ "Virulent",
    Temperate >= threshold ~ "Temperate",
    .default = "uncertain"
  ))

## PLOT

subset <- c("all", "core", "non-core")

for (set in subset) {
  filt_tbl <- bacphlip_assignment
  if (set == "core") {
    filt_tbl <- filt_tbl %>%
      filter(contig %in% present_in_all_countries)
  }
  if (set == "non-core") {
    filt_tbl <- filt_tbl %>%
      filter(!contig %in% present_in_all_countries)
  }
}

updated_classif <- classification %>%
  left_join(., bacphlip_assignment, by = "contig") %>%
  mutate(Lifestyle_bacphlip = assignment) %>% 
  relocate(Lifestyle_bacphlip, .after = Core) %>%
  select(-c(assignment, Virulent, Temperate))

bacphlip_lifestyle_plot <- updated_classif %>%
  filter(
    Class == "Caudoviricetes",
    Lifestyle_bacphlip != "uncertain"
  ) %>%
  count(Core, Lifestyle_bacphlip) %>%
  group_by(Core) %>%
  mutate(
    prop = n / sum(n),
    Core = factor(Core, levels = c("yes", "no")),
    Lifestyle_bacphlip = factor(Lifestyle_bacphlip, levels = c("Virulent", "Temperate"))
  ) %>%
  ungroup() %>%
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

## For convenience and to avoid back-tracking:
# write_csv(updated_classif, "output/R/classification.csv")

