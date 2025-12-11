library(tidyverse)
library(patchwork)
library(ggtext)

source("scripts/R/helpers/mixed_helpers.R")


#### Main figure
plots <- list()
plots$gene_presence$`PAPS reductase`$Insecticides <- readRDS("output/R/genes_pathogens_and_landuse/gene_presence_vs_landuse/single_panels/RDS.PAPS reductase.Insecticides.rds") +
  scale_y_continuous(limits = c(NA, NA),
                     breaks = c(0.76, 0.92),
                     labels = c(0.76, 0.92))
plots$pathogen_ct$`DWV B`$Cropland_in_2km_radius <- readRDS("output/R/genes_pathogens_and_landuse/pathogen_ct_vs_landuse/single_panels/RDS.DWV B.Cropland_in_2km_radius.rds") +
  scale_y_continuous(limits = c(NA, NA),
                     breaks = c(24.84, 33.22),
                     labels = c(24.84, 33.22))
plots$pathogen_ct$BQCV$Insecticides <- readRDS("output/R/genes_pathogens_and_landuse/pathogen_ct_vs_landuse/single_panels/RDS.BQCV.Insecticides.rds") +
  scale_y_continuous(limits = c(NA, NA),
                     breaks = c(23.37, 26.53),
                     labels = c(23.37, 26.53))
plots$gene_tpm$`PAPS reductase`$Cropland_in_2km_radius <- readRDS("output/R/genes_pathogens_and_landuse/gene_tpm_vs_landuse/single_panels/RDS.PAPS reductase.Cropland_in_2km_radius.rds")
plots$gene_tpm$`PAPS reductase`$`Pesticides (total)` <- readRDS("output/R/genes_pathogens_and_landuse/gene_tpm_vs_landuse/single_panels/RDS.PAPS reductase.Pesticides (total).rds")
plots$gene_tpm$`PAPS reductase`$`Fungicides and Bactericides` <- readRDS("output/R/genes_pathogens_and_landuse/gene_tpm_vs_landuse/single_panels/RDS.PAPS reductase.Fungicides and Bactericides.rds")
plots$gene_tpm$`PAPS reductase`$Herbicides <- readRDS("output/R/genes_pathogens_and_landuse/gene_tpm_vs_landuse/single_panels/RDS.PAPS reductase.Herbicides.rds")
plots$gene_tpm$`PAPS reductase`$`Plant Growth Regulators` <- readRDS("output/R/genes_pathogens_and_landuse/gene_tpm_vs_landuse/single_panels/RDS.PAPS reductase.Plant Growth Regulators.rds")

color_list <- list(dark = list(
  `PAPS reductase` = "#ef8f01",
  `DWV B` = "#6A0DAD",
  BQCV = "black",
  SBV = "#2A6666"
  ),
  bright = list(
    `PAPS reductase` = "#ef8f01",
    `DWV B` = "#6A0DAD",
    BQCV = "black",
    SBV = "#2A6666"
    )
  )

common_legend_1 <- legend_factory(title = "Gene", 
                                  items = names(color_list$dark)[1],
                                  colors = unlist(color_list$dark)[1],
                                  position = "top")
common_legend_1b <- legend_factory(title = "Gene", 
                                  items = names(color_list$dark)[1],
                                  colors = unlist(color_list$dark)[1],
                                  position = "left")
common_legend_2 <- legend_factory(title = "Pathogen",
                                     items = names(color_list$dark)[c(2,3)],
                                     colors = unlist(color_list$dark)[c(2,3)],
                                     position = "top")
common_legend_2b <- legend_factory(title = "Pathogen",
                                  items = names(color_list$dark)[c(2,3)],
                                  colors = unlist(color_list$dark)[c(2,3)],
                                  position = "left")
common_legend_3 <- legend_factory(title = "Pathogen",
                                   items = names(color_list$dark)[c(2:4)],
                                   colors = unlist(color_list$dark)[c(2:4)],
                                   position = "left")
# common_legend_row3b <- legend_factory(title = "Gene",
#                                       items = names(color_list$dark)[1],
#                                       colors = unlist(color_list$dark)[1],
#                                       position = "left") /
#   plot_spacer() + plot_layout(heights = c(25,1))

paps_square <- plots$gene_tpm$`PAPS reductase`$Cropland_in_2km_radius +
  plots$gene_tpm$`PAPS reductase`$`Pesticides (total)` +
  plots$gene_tpm$`PAPS reductase`$`Fungicides and Bactericides` +
  plots$gene_tpm$`PAPS reductase`$Herbicides +
  plot_spacer() +
  plots$gene_tpm$`PAPS reductase`$`Plant Growth Regulators` +
  plot_layout(ncol = 2, axes = "collect")

pathogens <-   
  plots$pathogen_ct$`DWV B`$Cropland_in_2km_radius +
  plots$pathogen_ct$BQCV$Insecticides +
  plot_layout(nrow = 1, axes = "collect")

#### Supplementary figure

names_presence <- c(
  "RDS.PAPS reductase.Herbicides – Dinitroanilines.rds",
  "RDS.PAPS reductase.Insecticides - nes.rds"
)

names_tpm <- c(
  "RDS.PAPS reductase.Fung & Bact – Dithiocarbamates.rds",
  "RDS.PAPS reductase.Fung & Bact - nes.rds",
  "RDS.PAPS reductase.Herbicides – Amides.rds",
  "RDS.PAPS reductase.Herbicides – Phenoxy hormone products.rds",
  "RDS.PAPS reductase.Herbicides – Triazines.rds",
  "RDS.PAPS reductase.Herbicides – Urea derivates.rds",
  "RDS.PAPS reductase.Herbicides - nes.rds",
  "RDS.PAPS reductase.Insecticides – Carbamates.rds",
  "RDS.PAPS reductase.Insecticides – Pyrethroids.rds",
  "RDS.PAPS reductase.Fung & Bact – Benzimidazoles.rds",
  "RDS.PAPS reductase.Fung & Bact – Diazines, morpholines.rds",
  "RDS.PAPS reductase.Fung & Bact – Triazoles, diazoles.rds",
  "RDS.PAPS reductase.Herbicides – Carbamates.rds"
)
names_ct <- c(
  "RDS.BQCV.Insecticides - nes.rds",
  "RDS.DWV B.Herbicides – Dinitroanilines.rds",
  "RDS.SBV.Herbicides – Phenoxy hormone products.rds"
)


supp_plots <- list()

for (file in names_presence) {
  supp_plots$presence[[file]] <- readRDS(paste0("output/R/genes_pathogens_and_landuse/gene_presence_vs_landuse/single_panels/", file))
}
for (file in names_tpm) {
  supp_plots$tpm[[file]] <- readRDS(paste0("output/R/genes_pathogens_and_landuse/gene_tpm_vs_landuse/single_panels/", file))
}
for (file in names_ct) {
  supp_plots$ct[[file]] <- readRDS(paste0("output/R/genes_pathogens_and_landuse/pathogen_ct_vs_landuse/single_panels/", file))
}

sup_figure_presence <- wrap_plots(supp_plots$presence, axes = "collect_y")
sup_figure_tpm <- wrap_plots(supp_plots$tpm, axes = "collect_y")
sup_figure_ct <- wrap_plots(supp_plots$ct, axes = "collect_y")

#####
# Save files
system("mkdir -p output/R/genes_pathogens_and_landuse/selected_graphs")

# Size of curve plots: per curve high 2.5, wide 2.5 then total width +1

# Main figure
ggsave("output/R/genes_pathogens_and_landuse/selected_graphs/main_paps_tpm.pdf",
       paps_square, height = 8.5, width = 6)
ggsave("output/R/genes_pathogens_and_landuse/selected_graphs/main_paps_presence.pdf",
       plots$gene_presence$`PAPS reductase`$Insecticides, height = 2.83, width = 3)
ggsave("output/R/genes_pathogens_and_landuse/selected_graphs/main_pathogen_ct.pdf",
       pathogens, height = 2.83, width = 6)


# Legends
ggsave("output/R/genes_pathogens_and_landuse/selected_graphs/legend_paps_top.pdf",
       common_legend_1, height = 1, width = 3)
ggsave("output/R/genes_pathogens_and_landuse/selected_graphs/legend_paps_left.pdf",
       common_legend_1b, height = 1, width = 3)
ggsave("output/R/genes_pathogens_and_landuse/selected_graphs/legend_viruses_top.pdf",
       common_legend_2, height = 1, width = 3)
ggsave("output/R/genes_pathogens_and_landuse/selected_graphs/legend_viruses_left.pdf",
       common_legend_2b, height = 1, width = 3)
ggsave("output/R/genes_pathogens_and_landuse/selected_graphs/legend_viruses_left_all.pdf",
       common_legend_3, height = 1, width = 3)

# Supp figure
ggsave("output/R/genes_pathogens_and_landuse/selected_graphs/supp_paps_presence.pdf",
       sup_figure_presence, height = 3, width = 6)
ggsave("output/R/genes_pathogens_and_landuse/selected_graphs/supp_paps_tpm.pdf",
       sup_figure_tpm, height = 12, width = 12)
ggsave("output/R/genes_pathogens_and_landuse/selected_graphs/supp_pathogen_ct.pdf",
       sup_figure_ct, height = 3, width = 9)

